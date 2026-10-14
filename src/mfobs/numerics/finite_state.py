"""Finite-state projection of an exported constant-jump event process.

No graph rules, pattern names or target moment equations are supplied here.
The query must be a nondecreasing count absent from every propensity. We sum
out that coordinate exactly, carrying its first two *unnormalised* conditional
moments over the remaining states. With forward matrix A and reward matrices
B_k[i,j] = sum_{e:j->i} rate_e(j) * delta_query(e)**k, the equations are

    p' = A p,   m' = A m + B_1 p,   s' = A s + 2 B_1 m + B_2 p.

Outgoing transitions kill probability; they are never reflected or renormalised.
These are finite-state forward equations, not a finite moment-closure claim.

Supported truncation certificate: non-Boolean state coordinates may increase
only through constant-rate events. A linear envelope q or q+z_i must have the
same property. These conditions are checked from the event expressions and
jumps. Independent constant-rate event clocks give Poisson upper bounds on
the coordinates and query. A union bound controls exit probability; Cauchy-
Schwarz with envelope moments of orders 2 and 4 controls the missing first and
second moments, respectively. No independence between exit and query is assumed.
Unsupported models fail explicitly rather than receiving an uncertified bound.

The reported bounds concern state truncation. Sparse matrix exponentiation is
performed in floating-point arithmetic; its rounding error is not certified.
FSP: Munsky & Khammash (2006), doi:10.1063/1.2145882.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from math import ceil
from typing import Mapping, Sequence

import numpy as np
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.stats import poisson
import sympy as sp


@dataclass
class MomentSolution:
    time: np.ndarray
    mean: np.ndarray
    variance: np.ndarray
    second_moment: np.ndarray
    surviving_probability: np.ndarray
    diagnostics: dict


def _poisson_envelope_moments(offset, step, intensity):
    """Raw moments of offset + step * Poisson(intensity), through order four."""
    lam = float(intensity)
    raw = [1., lam, lam**2 + lam, lam**3 + 3*lam**2 + lam,
           lam**4 + 6*lam**3 + 7*lam**2 + lam]
    from math import comb
    return [sum(comb(k, j)*offset**(k-j)*step**j*raw[j] for j in range(k+1))
            for k in range(5)]


def _certificate(names, deltas, expressions, initial, query_index, booleans,
                 horizon, tolerance):
    """Prove immigration envelopes by checking every event jump symbolically."""
    constant = np.array([not a.free_symbols for a in expressions], dtype=bool)
    constant_rates = np.array([float(a) if c else 0.
                               for a, c in zip(expressions, constant)])
    if np.any(constant_rates < 0) or not np.all(np.isfinite(constant_rates)):
        raise ValueError('Constant rates must be finite and nonnegative')
    # Zero-rate events cannot affect a bound or the reachable states.
    active = np.array([a != 0 for a in expressions], dtype=bool)
    resources = [i for i, name in enumerate(names)
                 if i != query_index and name not in booleans]
    immigration = {}
    for i in resources:
        positive = (deltas[:, i] > 0) & active
        if np.any(positive & ~constant):
            raise ValueError(f'No immigration truncation certificate for {names[i]}: '
                             'a state-dependent event increases it')
        immigration[i] = (int(max(deltas[positive, i], default=1)),
                          float(constant_rates[positive].sum()) * horizon)
    # q <= q + z_i because all coordinates are nonnegative. Seek a witness,
    # rather than using a pool name or the selected word's final letter.
    candidates = []
    for extra in [None, *resources]:
        weights = np.zeros(len(names), dtype=np.int64)
        weights[query_index] = 1
        if extra is not None:
            weights[extra] = 1
        jump = deltas @ weights
        positive = (jump > 0) & active
        if np.any(positive & ~constant):
            continue
        step = int(max(jump[positive], default=1))
        intensity = float(constant_rates[positive].sum()) * horizon
        offset = int(initial @ weights)
        moments = _poisson_envelope_moments(offset, step, intensity)
        candidates.append((moments[4], weights, offset, step, intensity, moments))
    if not candidates:
        raise ValueError('No query-moment envelope found among q and q + one coordinate')
    _, weights, offset, step, intensity, moments = min(candidates, key=lambda x: x[0])

    def errors(probability):
        e1 = float(np.sqrt(moments[2] * probability))
        e2 = float(np.sqrt(moments[4] * probability))
        # |(s2-mu^2)-(s2_trunc-mu_trunc^2)| <= e2 + 2 E[B(T)] e1.
        return e1, e2, e2 + 2*moments[1]*e1

    # Avoid inverse survival functions, which lose tiny tail probabilities.
    counts = {i: int(ceil(lam)) for i, (_, lam) in immigration.items()}
    while True:
        tails = {i: float(poisson.sf(counts[i], lam))
                 for i, (_, lam) in immigration.items()}
        probability = min(1., sum(tails.values()))
        e1, e2, ev = errors(probability)
        if max(e1, ev) <= tolerance:
            break
        if not tails or max(counts.values(), default=0) > 10000:
            raise RuntimeError('Truncation bounds could not be met within the cutoff budget')
        counts[max(tails, key=tails.get)] += 1
    bounds = {names[i]: int(initial[i] + immigration[i][0]*n)
              for i, n in counts.items()}
    bounds.update({name: 1 for name in booleans})
    envelope = {name: int(w) for name, w in zip(names, weights) if w}
    return bounds, dict(
        cutoff=bounds, horizon=float(horizon), truncation_tolerance=float(tolerance),
        exit_probability_bound=probability, mean_error_bound=e1,
        second_moment_error_bound=e2, variance_error_bound=ev,
        query_envelope=dict(weights=envelope, initial=offset, poisson_step=step,
                            poisson_mean_at_horizon=intensity),
        immigration={names[i]: {'step': jump, 'poisson_mean_at_horizon': lam}
                     for i, (jump, lam) in immigration.items()},
        bound_method='Poisson immigration envelopes, union bound and Cauchy-Schwarz',
        floating_point_error='Not included in the analytical truncation bounds',
    )


def solve_moments(*, names: Sequence[str], deltas, propensities: Sequence[sp.Expr],
                  initial: Mapping[str, int], query: str, boolean_names: Sequence[str],
                  times, truncation_tolerance=1e-6, max_states=150000,
                  adaptive=None) -> MomentSolution:
    """Compute a passive count's mean/variance from numerical event propensities.

    Input propensities already include numerical rate constants. All counts are
    nonnegative integers; Boolean coordinates come from the certified model.
    The time grid must be uniform, start at zero and contain at least two points.
    The result is not conditioned on staying within the truncated state space.
    For large a-priori boxes, adaptive FSP integrates boundary outflow and its
    query envelope moments. It enlarges the box until their bounds meet the same
    tolerance. These bounds also exclude floating-point integration error.
    """
    names = tuple(names)
    symbols = tuple(sp.Symbol(n) for n in names)
    expressions = tuple(sp.expand(a) for a in propensities)
    raw_delta = np.asarray(deltas)
    time = np.asarray(times, dtype=float)
    if (not names or len(set(names)) != len(names) or query not in names
            or set(initial) != set(names) or set(boolean_names) - set(names)):
        raise ValueError('Inconsistent coordinates, initial state or query')
    if (raw_delta.shape != (len(expressions), len(names))
            or not np.all(np.isfinite(raw_delta)) or np.any(raw_delta != np.floor(raw_delta))):
        raise ValueError('Supply one integer jump per event and coordinate')
    jumps = raw_delta.astype(np.int64)
    if any(not isinstance(v, (int, np.integer)) or isinstance(v, (bool, np.bool_))
           or v < 0 for v in initial.values()):
        raise ValueError('Initial counts must be nonnegative integers')
    if (time.ndim != 1 or len(time) < 2 or not np.all(np.isfinite(time))
            or time[0] != 0 or np.any(np.diff(time) <= 0)
            or not np.allclose(time, np.linspace(0, time[-1], len(time)), rtol=1e-12, atol=1e-14)):
        raise ValueError('Use a finite, increasing, uniform time grid starting at zero')
    if not np.isfinite(truncation_tolerance) or truncation_tolerance <= 0 or max_states < 1:
        raise ValueError('Use a positive tolerance and state budget')
    if any(a.free_symbols - set(symbols) for a in expressions):
        raise ValueError('Numerical propensities still contain unknown parameters')
    qi = names.index(query)
    if any(symbols[qi] in a.free_symbols for a in expressions) or np.any(jumps[:, qi] < 0):
        raise ValueError('The query must be a nondecreasing count absent from all propensities')
    if query in boolean_names:
        raise ValueError('Use a count query, not a Boolean state coordinate')
    start = np.array([initial[n] for n in names], dtype=np.int64)
    if any(initial[n] not in (0, 1) for n in boolean_names):
        raise ValueError('Initial Boolean values must be 0 or 1')
    bounds, diagnostics = _certificate(names, jumps, expressions, start, qi,
                                       set(boolean_names), time[-1], truncation_tolerance)
    prior_bounds = dict(bounds)
    if adaptive is None:
        from math import prod
        adaptive = prod(v+1 for v in bounds.values()) > max_states
    if adaptive:
        for name, clock in diagnostics['immigration'].items():
            guess = initial[name] + clock['step']*int(ceil(8+8*np.sqrt(clock['poisson_mean_at_horizon'])))
            bounds[name] = min(bounds[name], guess)
    attempts = []
    while True:
        solution, exits = _project(names, symbols, jumps, expressions, start, qi,
                                   set(boolean_names), time, bounds, diagnostics, max_states)
        info = solution.diagnostics
        attempts.append(dict(cutoff=dict(bounds), states=info['states'],
                             mean_error_bound=info['mean_error_bound'],
                             variance_error_bound=info['variance_error_bound']))
        # Both estimates are valid; use the a-priori certificate when using its box.
        if bounds == prior_bounds:
            for field in ('mean_error_bound', 'second_moment_error_bound',
                          'variance_error_bound', 'exit_probability_bound'):
                info[field] = min(info[field], diagnostics[field])
            info['bound_method'] += '; capped by the a-priori Poisson bound'
        if max(info['mean_error_bound'], info['variance_error_bound']) <= truncation_tolerance:
            info.update(projection_attempts=attempts, max_states=max_states)
            return solution
        if not adaptive:
            raise RuntimeError('The chosen FSP did not meet the requested moment tolerance')
        available = [name for name in diagnostics['immigration'] if bounds[name] < prior_bounds[name]]
        if not available:
            raise RuntimeError('FSP moment error bounds did not converge')
        # Grow the face carrying most outgoing probability, rather than both axes.
        name = max(available, key=lambda n: exits.get(n, 0.))
        bounds[name] = min(prior_bounds[name], max(bounds[name]+8, int(ceil(1.35*bounds[name]))))


def _project(names, symbols, jumps, expressions, start, qi, booleans, time,
             bounds, certificate, max_states):
    """Forward equations plus first-exit envelope moments, with no reflection.

    At first exit, set U=q_after+w_base*z_after. Thereafter let U grow by the
    constant-rate Poisson envelope. Its probability-weighted first/second moments
    upper-bound the query moments on exited paths, including subsequent re-entry.
    The strong Markov property of the constant-rate clocks justifies the future
    immigration terms. Accumulating outflow directly avoids subtracting 1-mass.
    """
    base_indices = [i for i in range(len(names)) if i != qi]
    base_names = [names[i] for i in base_indices]
    base_delta = jumps[:, base_indices]
    caps = np.array([bounds[n] for n in base_names], dtype=np.int64)
    base_symbols = [symbols[i] for i in base_indices]
    evaluate = sp.lambdify(base_symbols, expressions, modules='numpy', cse=True)
    initial_base = tuple(map(int, start[base_indices]))
    states, indices = [initial_base], {initial_base: 0}
    pending = deque([0])
    row, col, values, r1, r2 = [], [], [], [], []
    envelope = certificate['query_envelope']
    weights = np.array([envelope['weights'].get(name, 0) for name in base_names])
    resources = list(certificate['immigration'])
    exit_col, exit_zero, exit_first, exit_second = [], [], [], []
    face_rows, face_cols, face_values = [], [], []
    while pending:
        j = pending.popleft()
        state = np.array(states[j], dtype=np.int64)
        rates = np.asarray(evaluate(*state), dtype=float).reshape(-1)
        if not np.all(np.isfinite(rates)) or np.any(rates < -1e-10):
            raise ValueError(f'Invalid event rate at {dict(zip(base_names, state))}')
        rates = np.maximum(rates, 0.)
        row.append(j); col.append(j); values.append(-float(rates.sum()))
        r1.append(0.); r2.append(0.)
        for e in np.flatnonzero(rates):
            target = state + base_delta[e]
            if np.any(target < 0) or any(target[i] > 1 for i, n in enumerate(base_names)
                                         if n in booleans):
                raise ValueError('An enabled event leaves the declared coordinate domain')
            if np.any(target > caps):
                offset = int(jumps[e, qi] + weights @ target)
                exit_col.append(j); exit_zero.append(rates[e])
                exit_first.append(rates[e]*offset); exit_second.append(rates[e]*offset**2)
                for i in np.flatnonzero(target > caps):
                    face_rows.append(resources.index(base_names[i]))
                    face_cols.append(j); face_values.append(rates[e])
                continue  # outflow remains on the diagonal: probability is killed
            key = tuple(map(int, target))
            if key not in indices:
                if len(states) >= max_states:
                    raise RuntimeError(f'FSP state budget {max_states:,} exhausted at t={time[-1]:g}; '
                                       'increase max_states or loosen truncation_tolerance. No curve returned.')
                indices[key] = len(states)
                pending.append(len(states))
                states.append(key)
            row.append(indices[key]); col.append(j); values.append(float(rates[e]))
            r1.append(float(rates[e] * jumps[e, qi]))
            r2.append(float(rates[e] * jumps[e, qi]**2))
    n = len(states)
    def matrix(entries):
        a = coo_matrix((entries, (row, col)), shape=(n, n)).tocsr()
        a.eliminate_zeros()
        return a
    A, B1, B2 = matrix(values), matrix(r1), matrix(r2)
    # Sink variables: exit mass, first/second envelope moments, face exit masses.
    k = 3+len(resources)
    def sink_column(entries, sink_row):
        return coo_matrix((entries, ([sink_row]*len(exit_col), exit_col)), shape=(k, n)).tocsr()
    from_p = (sink_column(exit_zero, 0) + sink_column(exit_first, 1)
              + sink_column(exit_second, 2)
              + coo_matrix((face_values, (np.array(face_rows, dtype=int)+3, face_cols)), shape=(k, n)).tocsr())
    from_m = sink_column(exit_zero, 1) + 2*sink_column(exit_first, 2)
    from_s = sink_column(exit_zero, 2)
    rate = envelope['poisson_mean_at_horizon']/time[-1]
    g, h = envelope['poisson_step']*rate, envelope['poisson_step']**2*rate
    sink = coo_matrix(([g, h, 2*g], ([1, 2, 2], [0, 0, 1])), shape=(k, k)).tocsr()
    operator = bmat([[A, None, None, None], [B1, A, None, None],
                     [B2, 2*B1, A, None], [from_p, from_m, from_s, sink]], format='csr')
    y0 = np.zeros(3*n+k)
    y0[0], y0[n], y0[2*n] = 1., float(start[qi]), float(start[qi]**2)
    # A diagonal similarity transform keeps the moment/sink blocks at comparable
    # scales. It changes no equations and avoids excessive exponential scaling.
    unit = max(1., envelope['initial'] + envelope['poisson_step']*envelope['poisson_mean_at_horizon'])
    scale = np.ones(3*n+k)
    scale[n:2*n] = unit; scale[2*n:3*n] = unit**2
    scale[3*n+1:3*n+3] = [unit, unit**2]
    operator = operator.multiply(scale[None, :]).multiply(1/scale[:, None]).tocsr()
    y0 /= scale
    # Stream a few time points per block, rather than retaining every full vector.
    totals = np.empty((len(time), 3+k))
    totals[0] = [1., float(start[qi]), float(start[qi]**2), *([0.]*k)]
    for first in range(1, len(time), 32):
        last = min(first+31, len(time)-1)
        y = expm_multiply(operator, y0, start=0., stop=float(time[last]-time[first-1]),
                          num=last-first+2, traceA=float(operator.diagonal().sum()))
        if not np.all(np.isfinite(y)) or np.min(y) < -1e-9:
            raise RuntimeError('Numerical solution failed positivity checks')
        totals[first:last+1, :3] = np.column_stack([y[1:, j*n:(j+1)*n].sum(axis=1)*unit**j for j in range(3)])
        totals[first:last+1, 3:] = y[1:, 3*n:]*scale[3*n:]
        y0 = y[-1].copy()
    mass, mean, second = totals[:, :3].T
    escaped, e1, e2 = totals[-1, 3:6]
    evar = e2 + 2*(envelope['initial'] + envelope['poisson_step']*envelope['poisson_mean_at_horizon'])*e1
    variance = second - mean**2
    if (np.min(variance) < -1e-8
            or np.max(mass) > 1+1e-9
            or np.max(np.abs(mass+totals[:, 3]-1)) > 1e-9):
        raise RuntimeError('Numerical solution failed positivity or probability-mass checks')
    diagnostics = {**certificate, 'cutoff': dict(bounds), 'mean_error_bound': float(e1),
                   'second_moment_error_bound': float(e2), 'variance_error_bound': float(evar),
                   'exit_probability_bound': float(escaped),
                   'bound_method': 'Integrated first-exit flux and Poisson query-envelope moments'}
    diagnostics.update(states=n, equations=3*n+k, query=names[qi],
                       base_coordinates=base_names, observed_mass_defect=float(max(0., 1-mass[-1])),
                       method='Finite-state projection; sparse matrix exponential; conditional query moments')
    return (MomentSolution(time, mean, np.maximum(variance, 0.), second, mass, diagnostics),
            dict(zip(resources, totals[-1, 6:])))
