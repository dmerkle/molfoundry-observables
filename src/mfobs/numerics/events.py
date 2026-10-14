"""Numerical consumers of exported event records, shared by scripts and notebooks."""
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping
import numpy as np
import pandas as pd
import sympy as sp

def analysis(record):
    return record.get("analysis", record)

def polynomial(record: dict, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    """Read sparse rational terms, never execute/sympify the stored text field."""
    result = sp.Integer(0)
    for term in record["terms"]:
        powers = term.get("monomial", term.get("powers"))
        if not isinstance(powers, dict):
            raise ValueError("Expected named monomial exponents")
        f = Fraction(term["coefficient"])
        value = sp.Rational(f.numerator, f.denominator)
        for name, power in powers.items():
            if name not in symbols or type(power) is not int or power < 0:
                raise ValueError(f"Unsupported polynomial factor {name}^{power}")
            value *= symbols[name] ** power
        result += value
    return sp.expand(result)


@dataclass
class EventSystem:
    names: tuple[str, ...]
    symbols: tuple[sp.Symbol, ...]
    identifiers: tuple[str, ...]
    channels: tuple[str, ...]
    deltas: np.ndarray
    activities: tuple[sp.Expr, ...]

    def drift(self) -> list[sp.Expr]:
        return [sp.expand(sum(sp.Symbol(rate) * delta[i] * activity for rate, delta, activity in zip(self.channels, self.deltas, self.activities))) for i in range(len(self.names))]


def event_system(record: dict) -> EventSystem:
    result = analysis(record)
    pjc = result.get("pjc", {})
    if pjc.get("status") != "CLOSED":
        raise ValueError("No certified finite event model in this record")
    names = tuple(pjc["observable_coordinates"])
    if not names or len(names) != len(set(names)):
        raise ValueError("Invalid coordinate list")
    symbols = {name: sp.Symbol(name) for name in names}
    rate_symbols = {r["id"]: r["rate"] for r in result.get("contract", {}).get("normalized_rules", [])}
    terms = pjc["operator"]["terms"]
    rows, acts = [], []
    for term in terms:
        delta = term["delta"]
        if set(delta) - set(names) or any(type(v) is not int for v in delta.values()):
            raise ValueError("Event has an unknown coordinate or a noninteger jump")
        if rate_symbols.get(term["rate"], term["rate"]) in symbols:
            raise ValueError("Rate and coordinate names collide")
        rows.append([delta.get(name, 0) for name in names])
        acts.append(polynomial(term["propensity_polynomial"], symbols))
    return EventSystem(names, tuple(symbols.values()), tuple(t["event_class"] for t in terms), tuple(rate_symbols.get(t["rate"], t["rate"]) for t in terms), np.array(rows, dtype=np.int64).reshape(-1, len(names)), tuple(acts))


def numeric_parameters(events: EventSystem, parameters: Mapping[str, float]) -> dict:
    required = set(events.channels)
    if set(parameters) != required:
        raise ValueError(f"Supply exactly these rate parameters: {sorted(required)}; got {sorted(parameters)}")
    if any(not np.isfinite(v) or v < 0 for v in parameters.values()):
        raise ValueError("Rate parameters must be finite and nonnegative")
    return {sp.Symbol(k): float(v) for k, v in parameters.items()}


def affine_generator(record: dict, parameters: Mapping[str, float]):
    """Build the affine drift from recorded events; refuse unclosed moments."""
    events = event_system(record)
    subs = numeric_parameters(events, parameters)
    zero = {s: 0 for s in events.symbols}
    A, b = [], []
    for expr in events.drift():
        poly = sp.Poly(expr, *events.symbols)
        if poly.total_degree() > 1:
            raise ValueError("Coordinate drift is nonlinear. Its first-moment identities cannot be integrated as a closed expectation ODE.")
        A.append([float(sp.diff(expr, s).subs(subs)) for s in events.symbols])
        b.append(float(expr.subs(zero).subs(subs)))
    return events.names, np.asarray(A), np.asarray(b)


def integrate_events(record: dict, initial: Mapping[str, int], parameters: Mapping[str, float],
                     times, *, query: str, truncation_tolerance=1e-6, max_states=150000):
    """Numerical mean/variance of a count from the exported events, with FSP bounds."""
    from mfobs.numerics.finite_state import solve_moments
    events = event_system(record)
    subs = numeric_parameters(events, parameters)
    return solve_moments(
        names=events.names, deltas=events.deltas,
        propensities=[sp.expand(sp.Symbol(c).subs(subs)*a)
                      for c, a in zip(events.channels, events.activities)],
        initial=initial, query=query,
        boolean_names=analysis(record)['pjc'].get('quotient_identities', {}).get('boolean', []),
        times=times, truncation_tolerance=truncation_tolerance, max_states=max_states,
    )


def simulate_events(record: dict, initial: Mapping[str, int], parameters: Mapping[str, float], times, rng, *, max_events=100000):
    """Direct SSA on the recorded finite-event generator; no mean-field step.

    Negative/nonfinite rates, invalid counts and exhausted budgets raise errors.
    This numerical sampling does not itself verify the model's reduction theorem.
    The returned trace keeps the first 25 events; ``trace.attrs["event_count"]``
    gives the total number of events simulated.
    """
    events = event_system(record)
    subs = numeric_parameters(events, parameters)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not len(times) or times[0] != 0 or np.any(np.diff(times) <= 0) or not np.all(np.isfinite(times)):
        raise ValueError("times must be finite, strictly increasing, and start at zero")
    if set(initial) != set(events.names) or any(type(v) not in (int, np.int64, np.int32) or v < 0 for v in initial.values()):
        raise ValueError("Supply a nonnegative integer initial value for every coordinate")
    if max_events < 1:
        raise ValueError("max_events must be positive")
    evaluate = sp.lambdify(events.symbols, [sp.Symbol(c).subs(subs) * a for c, a in zip(events.channels, events.activities)], modules="numpy")
    state = np.array([initial[n] for n in events.names], dtype=np.int64)
    samples = np.empty((len(times), len(state)), dtype=np.int64)
    clock, cursor, fired = 0., 0, 0
    trace = []
    while cursor < len(times):
        rates = np.asarray(evaluate(*state), dtype=float).reshape(-1)
        if len(rates) != len(events.identifiers) or not np.all(np.isfinite(rates)) or np.any(rates < -1e-10):
            raise ValueError(f"Invalid rates at {dict(zip(events.names, state))}: {rates}")
        rates = np.maximum(rates, 0.)  # only tolerate floating-point roundoff
        total = rates.sum()
        if not np.isfinite(total):
            raise ValueError("Total event rate is not finite")
        next_time = clock + rng.exponential(1 / total) if total else np.inf
        while cursor < len(times) and times[cursor] < next_time:
            samples[cursor] = state
            cursor += 1
        if cursor == len(times):
            break
        if fired >= max_events:
            raise RuntimeError("SSA event budget exhausted; no truncated trajectory returned")
        which = int(rng.choice(len(rates), p=rates / total))
        state = state + events.deltas[which]
        if np.any(state < 0):
            raise ValueError("An enabled jump produced a negative pattern count")
        clock = next_time
        fired += 1
        if len(trace) < 25:
            trace.append({"time": clock, "event": events.identifiers[which], **dict(zip(events.names, map(int, state)))})
    trace = pd.DataFrame(trace)
    trace.attrs["event_count"] = fired
    return samples, trace



def sample_moments(record, initial, parameters, times, *, paths=100, seed=20260924,
                   truncation_tolerance=1e-6, max_states=150000):
    """FSP moments and independent SSA samples from the same exported events."""
    if not isinstance(paths, int) or paths < 2:
        raise ValueError('Use at least two SSA paths')
    query = analysis(record)['support_discovery']['selected_queries'][0]
    solution = integrate_events(record, initial, parameters, times, query=query,
        truncation_tolerance=truncation_tolerance, max_states=max_states)
    seeds = np.random.SeedSequence(seed).spawn(paths)
    samples = np.stack([simulate_events(record, initial, parameters, times,
                       np.random.default_rng(s))[0] for s in seeds])
    column = event_system(record).names.index(query)
    values = samples[:, :, column]
    return dict(time=np.asarray(times), samples=samples, mean=solution.mean,
        variance=solution.variance, second_moment=solution.second_moment,
        ssa_mean=values.mean(0), ssa_variance=values.var(0, ddof=1),
        ssa_mean_error=1.96*values.std(0, ddof=1)/np.sqrt(paths),
        numerical=solution.diagnostics, rates=dict(parameters), ssa_paths=paths, ssa_seed=seed,
        options=dict(end=float(times[-1]), paths=paths, seed=seed,
                     truncation_tolerance=truncation_tolerance, max_states=max_states))
