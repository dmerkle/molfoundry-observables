"""Generic exact polynomial jump closure on a certified event quotient.

The model-specific input is finite: observable names and semantic event
classes.  Every event class supplies a constant observable jump and a
polynomial propensity.  The implementation then derives, without equations
hard-coded by an example,

* the rate-separated first-generator action;
* a FIFO linear-observable worklist;
* the finite span of all nested diagonal commutators; and
* the CTMC exponential-moment generating-function operator.

Finiteness of this table alone is not an all-state proof.  Callers must attach
an event-partition certificate obtained from their state language.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Iterable, Mapping, Sequence

from mfobs.core.operator_span import discover_diagonal_commutator_span


Monomial = tuple[int, ...]


def _fraction(value: int | str | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


@dataclass(frozen=True)
class Polynomial:
    """Canonical sparse polynomial in one fixed :class:`PolynomialRing`."""

    terms: tuple[tuple[Monomial, Fraction], ...]

    def as_dict(self) -> dict[Monomial, Fraction]:
        return dict(self.terms)

    @property
    def degree(self) -> int:
        return max((sum(monomial) for monomial, _ in self.terms), default=0)

    @property
    def is_zero(self) -> bool:
        return not self.terms


class PolynomialRing:
    """Sparse rational polynomials modulo declared Boolean state identities."""

    def __init__(
        self,
        variables: Sequence[str],
        *,
        boolean_variables: Iterable[str] = (),
        exclusive_groups: Sequence[Sequence[str]] = (),
        boolean_products: Mapping[tuple[str, str], str | None] | None = None,
    ) -> None:
        self.variables = tuple(variables)
        if not self.variables or len(self.variables) != len(set(self.variables)):
            raise ValueError("polynomial variables must be nonempty and unique")
        self.indices = {name: index for index, name in enumerate(self.variables)}
        try:
            self.boolean_indices = frozenset(
                self.indices[name] for name in boolean_variables
            )
            self.exclusive_groups = tuple(
                frozenset(self.indices[name] for name in group)
                for group in exclusive_groups
            )
        except KeyError as exc:
            raise ValueError(f"unknown quotient variable {exc.args[0]!r}") from exc
        if any(not group for group in self.exclusive_groups):
            raise ValueError("exclusive Boolean groups may not be empty")
        if any(not group <= self.boolean_indices for group in self.exclusive_groups):
            raise ValueError("exclusive groups must contain Boolean variables")
        self.boolean_products: dict[tuple[int, int], int | None] = {}
        for pair, target in (boolean_products or {}).items():
            if len(pair) != 2 or pair[0] not in self.indices or pair[1] not in self.indices:
                raise ValueError(f"invalid Boolean product relation {pair!r}")
            indices = tuple(sorted((self.indices[pair[0]], self.indices[pair[1]])))
            if not set(indices) <= self.boolean_indices:
                raise ValueError("Boolean product relation uses a non-Boolean variable")
            if target is not None:
                if target not in self.indices or self.indices[target] not in self.boolean_indices:
                    raise ValueError("Boolean product target must be Boolean")
                target_index = self.indices[target]
            else:
                target_index = None
            self.boolean_products[indices] = target_index

    @property
    def zero_monomial(self) -> Monomial:
        return (0,) * len(self.variables)

    def _reduce_monomial(self, monomial: Monomial) -> Monomial | None:
        if len(monomial) != len(self.variables) or any(
            type(exponent) is not int or exponent < 0 for exponent in monomial
        ):
            raise ValueError("invalid monomial exponent vector")
        reduced = list(monomial)
        for index in self.boolean_indices:
            reduced[index] = int(bool(reduced[index]))
        active = {index for index in self.boolean_indices if reduced[index]}
        changed = True
        while changed:
            changed = False
            for left in sorted(active):
                for right in sorted(active):
                    if left >= right:
                        continue
                    pair = (left, right)
                    if pair not in self.boolean_products:
                        continue
                    target = self.boolean_products[pair]
                    if target is None:
                        return None
                    active.remove(left)
                    active.remove(right)
                    active.add(target)
                    changed = True
                    break
                if changed:
                    break
        for index in self.boolean_indices:
            reduced[index] = int(index in active)
        for group in self.exclusive_groups:
            if sum(bool(reduced[index]) for index in group) > 1:
                return None
        return tuple(reduced)

    def polynomial(
        self, terms: Mapping[Monomial, int | str | Fraction] | None = None
    ) -> Polynomial:
        combined: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
        for monomial, coefficient in (terms or {}).items():
            reduced = self._reduce_monomial(tuple(monomial))
            if reduced is not None:
                combined[reduced] += _fraction(coefficient)
        return Polynomial(tuple(sorted(
            (monomial, coefficient)
            for monomial, coefficient in combined.items()
            if coefficient
        )))

    def constant(self, value: int | str | Fraction) -> Polynomial:
        return self.polynomial({self.zero_monomial: value})

    def coordinate(self, name: str) -> Polynomial:
        exponent = [0] * len(self.variables)
        exponent[self.indices[name]] = 1
        return self.polynomial({tuple(exponent): 1})

    def monomial(self, powers: Mapping[str, int]) -> Polynomial:
        exponent = [0] * len(self.variables)
        for name, power in powers.items():
            exponent[self.indices[name]] = power
        return self.polynomial({tuple(exponent): 1})

    def add(self, *values: Polynomial) -> Polynomial:
        terms: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
        for value in values:
            for monomial, coefficient in value.terms:
                terms[monomial] += coefficient
        return self.polynomial(terms)

    def scale(self, value: Polynomial, coefficient: int | str | Fraction) -> Polynomial:
        factor = _fraction(coefficient)
        return self.polynomial({
            monomial: factor * old for monomial, old in value.terms
        })

    def subtract(self, left: Polynomial, right: Polynomial) -> Polynomial:
        return self.add(left, self.scale(right, -1))

    def multiply(self, *values: Polynomial) -> Polynomial:
        result = self.constant(1)
        for value in values:
            terms: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
            for left, left_coefficient in result.terms:
                for right, right_coefficient in value.terms:
                    monomial = tuple(a + b for a, b in zip(left, right))
                    reduced = self._reduce_monomial(monomial)
                    if reduced is not None:
                        terms[reduced] += left_coefficient * right_coefficient
            result = self.polynomial(terms)
        return result

    def shift(self, value: Polynomial, delta: Sequence[int]) -> Polynomial:
        """Return ``value(x + delta)`` exactly by binomial expansion."""

        if len(delta) != len(self.variables):
            raise ValueError("jump dimension differs from polynomial ring")
        total = self.constant(0)
        for monomial, coefficient in value.terms:
            expanded = self.constant(coefficient)
            for index, exponent in enumerate(monomial):
                if not exponent:
                    continue
                terms = {}
                for retained_power in range(exponent + 1):
                    constant_power = exponent - retained_power
                    factor = comb(exponent, retained_power) * (
                        delta[index] ** constant_power
                    )
                    powers = [0] * len(self.variables)
                    powers[index] = retained_power
                    terms[tuple(powers)] = factor
                expanded = self.multiply(expanded, self.polynomial(terms))
            total = self.add(total, expanded)
        return total

    def evaluate(self, value: Polynomial, coordinates: Mapping[str, int]) -> Fraction:
        if set(coordinates) != set(self.variables):
            raise ValueError("evaluation must supply every coordinate exactly once")
        return sum((
            coefficient
            * self._evaluate_monomial(monomial, coordinates)
            for monomial, coefficient in value.terms
        ), Fraction())

    def _evaluate_monomial(
        self, monomial: Monomial, coordinates: Mapping[str, int]
    ) -> int:
        result = 1
        for name, exponent in zip(self.variables, monomial):
            result *= coordinates[name] ** exponent
        return result

    def monomial_text(self, monomial: Monomial) -> str:
        factors = []
        for name, exponent in zip(self.variables, monomial):
            if exponent == 1:
                factors.append(name)
            elif exponent:
                factors.append(f"{name}^{exponent}")
        return "*".join(factors) or "1"

    def text(self, value: Polynomial) -> str:
        pieces: list[str] = []
        for monomial, coefficient in value.terms:
            sign = "-" if coefficient < 0 else "+"
            magnitude = abs(coefficient)
            body = self.monomial_text(monomial)
            if body == "1":
                term = str(magnitude)
            elif magnitude == 1:
                term = body
            else:
                term = f"{magnitude}*{body}"
            if not pieces:
                pieces.append(term if sign == "+" else f"-{term}")
            else:
                pieces.append(f" {sign} {term}")
        return "".join(pieces) or "0"

    def record(self, value: Polynomial) -> dict:
        return {
            "text": self.text(value),
            "degree": value.degree,
            "terms": [
                {
                    "coefficient": str(coefficient),
                    "monomial": {
                        name: exponent
                        for name, exponent in zip(self.variables, monomial)
                        if exponent
                    },
                    "differential_monomial": self.monomial_text(monomial),
                }
                for monomial, coefficient in value.terms
            ],
        }


@dataclass(frozen=True)
class JumpEvent:
    """One all-state event class with a constant observable jump."""

    identifier: str
    channel: str
    delta: tuple[int, ...]
    propensity: Polynomial
    rule: str
    predicate: str


@dataclass(frozen=True)
class JumpModelProvenance:
    """Uniform origin record for every exact polynomial jump backend."""

    backend: str
    event_partition_origin: str
    jump_origin: str
    activity_origin: str
    jump_projection: str
    composition_scope: str
    all_jumps_composition_derived: bool
    supplied_delta_table: bool

    def __post_init__(self) -> None:
        text = (
            self.backend,
            self.event_partition_origin,
            self.jump_origin,
            self.activity_origin,
            self.jump_projection,
            self.composition_scope,
        )
        if any(not value.strip() for value in text):
            raise ValueError("jump-model provenance text must be nonempty")

    def record(self) -> dict:
        return {
            "schema_version": 1,
            "backend": self.backend,
            "event_partition_origin": self.event_partition_origin,
            "jump_origin": self.jump_origin,
            "activity_origin": self.activity_origin,
            "jump_projection": self.jump_projection,
            "composition_scope": self.composition_scope,
            "all_jumps_composition_derived": self.all_jumps_composition_derived,
            "supplied_delta_table": self.supplied_delta_table,
        }


@dataclass(frozen=True)
class PolynomialJumpModel:
    """Certified finite event quotient of an otherwise infinite state language."""

    name: str
    ring: PolynomialRing
    events: tuple[JumpEvent, ...]
    event_partition_certificate: Mapping
    input_definition_sha256: str
    provenance: JumpModelProvenance

    def __post_init__(self) -> None:
        identifiers = [event.identifier for event in self.events]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("event identifiers must be unique")
        for event in self.events:
            if len(event.delta) != len(self.ring.variables):
                raise ValueError(f"{event.identifier}: jump has the wrong dimension")
        if self.event_partition_certificate.get("status") != "VERIFIED_ALL_STATES":
            raise ValueError("pJC requires an all-state event-partition certificate")


def generator_action(
    model: PolynomialJumpModel,
    observable: Polynomial,
    *,
    channel: str | None = None,
) -> Polynomial:
    """Apply the exact backward generator to one polynomial observable."""

    ring = model.ring
    contributions = []
    for event in model.events:
        if channel is not None and event.channel != channel:
            continue
        jump = ring.subtract(ring.shift(observable, event.delta), observable)
        contributions.append(ring.multiply(event.propensity, jump))
    return ring.add(*contributions)


def one_step_event_expansion(
    model: PolynomialJumpModel,
    initial: Mapping[str, int],
) -> dict:
    """Expand one coordinate state through every enabled certified event class.

    This is a presentation/validation view of an already derived jump model,
    not an input to closure.  Activities and targets are evaluated directly
    from the generic polynomial ring and constant event jumps.
    """

    ring = model.ring
    if set(initial) != set(ring.variables):
        raise ValueError("initial state must supply every coordinate exactly once")
    source = {name: int(initial[name]) for name in ring.variables}
    transitions = []
    for event in model.events:
        activity = ring.evaluate(event.propensity, source)
        if activity <= 0:
            continue
        target = {
            name: source[name] + event.delta[index]
            for index, name in enumerate(ring.variables)
        }
        if any(value < 0 for value in target.values()):
            raise ValueError(
                f"{event.identifier}: positive activity leaves the nonnegative quotient"
            )
        transitions.append({
            "event_class": event.identifier,
            "rule": event.rule,
            "rate": event.channel,
            "delta": {
                name: value
                for name, value in zip(ring.variables, event.delta)
                if value
            },
            "activity": str(activity),
            "source": dict(source),
            "target": target,
        })
    return {
        "schema_version": 1,
        "source": source,
        "transitions": transitions,
        "event_partition_certificate": dict(model.event_partition_certificate),
        "input_definition_sha256": model.input_definition_sha256,
    }


def discover_linear_observables_from_events(
    model: PolynomialJumpModel,
    seeds: Sequence[str],
    *,
    maximum_dimension: int = 32,
) -> dict:
    """Linear closure for an explicitly supplied exact event generator."""
    from mfobs.core.linear_closure import close_linear_span

    ring = model.ring
    if not seeds or set(seeds) - set(ring.variables):
        raise ValueError("linear seeds must be polynomial-ring coordinates")
    expressions = {name: ring.coordinate(name) for name in dict.fromkeys(seeds)}

    def image(source, channel):
        value = generator_action(model, expressions[source], channel=channel)
        return value.as_dict(), {"generator_image": ring.record(value)}

    def designate(vector, index):
        # Names describe computed expressions; no preselected auxiliary basis.
        name = ring.text(ring.polynomial(vector))
        return name if name not in expressions else f"O{index}"

    def added(name, vector, source, channel, before, scale):
        expressions[name] = ring.polynomial(vector)

    result = close_linear_span(
        seeds={n: p.as_dict() for n, p in expressions.items()},
        channels=dict.fromkeys(event.channel for event in model.events),
        image=image, encode=lambda v: ring.record(ring.polynomial(v)),
        designate=designate, on_add=added, maximum_dimension=maximum_dimension,
    )
    return {
        "schema_version": 2, "status": result.status,
        "analysis": "linear-observable-discovery",
        "input_definition_sha256": model.input_definition_sha256,
        "basis": list(result.basis),
        "basis_expressions": {n: ring.record(p) for n, p in expressions.items()},
        "origins": result.origins, "steps": result.steps,
        "pending": result.pending, "obligations": result.obligations,
        "dimension": len(result.basis), "maximum_dimension": maximum_dimension,
        "scope": "exact on the certified state-language quotient",
    }



def _exponent_text(names: Sequence[str], delta: Sequence[int]) -> str:
    terms = []
    for name, jump in zip(names, delta):
        if not jump:
            continue
        variable = f"lambda_{name}"
        if jump == 1:
            terms.append(variable)
        elif jump == -1:
            terms.append(f"-{variable}")
        else:
            terms.append(f"{jump}*{variable}")
    return " + ".join(terms).replace("+ -", "- ") or "0"


def analyze_polynomial_jump_closure(model: PolynomialJumpModel) -> dict:
    """Derive the pJC certificate, nested-commutator span, and EMGF operator."""

    ring = model.ring
    jumps = {
        name: tuple(event.delta[index] for event in model.events)
        for index, name in enumerate(ring.variables)
    }
    span = discover_diagonal_commutator_span(
        event_classes=[event.identifier for event in model.events],
        jumps=jumps,
    )
    terms = []
    maximum_order = 0
    for event in model.events:
        maximum_order = max(maximum_order, event.propensity.degree)
        exponent = _exponent_text(ring.variables, event.delta)
        terms.append({
            "event_class": event.identifier,
            "rule": event.rule,
            "rate": event.channel,
            "predicate": event.predicate,
            "delta": {
                name: value for name, value in zip(ring.variables, event.delta) if value
            },
            "exponent": exponent,
            "factor": f"exp({exponent}) - 1",
            "propensity_polynomial": ring.record(event.propensity),
        })

    first_moments = {}
    for index, name in enumerate(ring.variables):
        by_channel = {}
        for channel in dict.fromkeys(event.channel for event in model.events):
            contributions = [
                ring.scale(event.propensity, event.delta[index])
                for event in model.events
                if event.channel == channel and event.delta[index]
            ]
            value = ring.add(*contributions)
            if not value.is_zero:
                by_channel[channel] = ring.record(value)
        first_moments[name] = by_channel

    return {
        "schema_version": 1,
        "status": "CLOSED",
        "analysis": "polynomial-jump-closure-emgf",
        "input_definition_sha256": model.input_definition_sha256,
        "scope": "all states in the certified state language; no finite-state inference",
        "event_partition_certificate": dict(model.event_partition_certificate),
        "jump_model_provenance": model.provenance.record(),
        "observable_coordinates": list(ring.variables),
        "quotient_identities": {
            "boolean": [ring.variables[index] for index in sorted(ring.boolean_indices)],
            "mutually_exclusive_groups": [
                [ring.variables[index] for index in sorted(group)]
                for group in ring.exclusive_groups
            ],
            "boolean_products": [
                {
                    "left": ring.variables[left],
                    "right": ring.variables[right],
                    "result": None if target is None else ring.variables[target],
                }
                for (left, right), target in sorted(ring.boolean_products.items())
            ],
        },
        "event_class_count": len(model.events),
        "operator_span": span.record(),
        "operator": {
            "formula": (
                "D = sum_e rate_e * (exp(lambda dot delta_e)-1) "
                "* a_e(partial_lambda)"
            ),
            "terms": terms,
            "maximum_differential_order": maximum_order,
        },
        "first_moment_generator_images": first_moments,
    }
