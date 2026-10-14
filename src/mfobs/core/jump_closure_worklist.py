"""Generic first-commutator and diagonal jump-projection synthesis.

The objects manipulated here are finite linear combinations of rules.  This
is important: an observable such as a typed-edge count is naturally a sum of
identity-rule components, not necessarily one graph.  Starting from supplied
identity expressions, the worklist repeatedly computes

    O  ->  [O, G]  ->  p([O, G]).

Here ``p`` is the DPO/SqPO jump projection of a first commutator; applying it
does not by itself assert polynomial jump closure of all nested commutators.
Any exact residual outside the current observable span becomes a new formal
identity expression and is queued.  The module knows nothing about trees,
chemistry, observable names, or published equations.  An application supplies
only a rule commutator (which may enforce its declared structural theory) and
an optional classifier used after discovery to give a result a readable name.
All arithmetic and span decisions are exact over the rationals.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction
from hashlib import sha256
from typing import Callable, Iterable, Mapping

from mfobs.core.engine import sqpo_jump_closure


def _fraction(value) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def _clean(values: Mapping[str, Fraction]) -> dict[str, Fraction]:
    return {key: _fraction(value) for key, value in values.items() if value}


def _encoded_coefficient(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _term_id(key: str) -> str:
    return sha256(key.encode()).hexdigest()


@dataclass(frozen=True)
class RuleExpression:
    """Sparse exact linear combination of canonical MolFoundry rules."""

    coefficients: Mapping[str, Fraction]
    representatives: Mapping[str, object]

    def __post_init__(self):
        coefficients = _clean(self.coefficients)
        if set(coefficients) - set(self.representatives):
            raise ValueError("every nonzero rule term needs a representative")
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(
            self,
            "representatives",
            {key: self.representatives[key] for key in coefficients},
        )

    @classmethod
    def from_rules(
        cls,
        rules: Iterable[object],
        coefficients: Iterable[int | Fraction] | None = None,
    ) -> "RuleExpression":
        rules = tuple(rules)
        values = tuple(coefficients) if coefficients is not None else (1,) * len(rules)
        if len(rules) != len(values):
            raise ValueError("one coefficient is required per rule")
        total: Counter = Counter()
        representatives = {}
        for rule, coefficient in zip(rules, values):
            key = rule.canonKey
            total[key] += _fraction(coefficient)
            representatives.setdefault(key, rule)
        return cls(total, representatives)

    def scale(self, coefficient: int | Fraction) -> "RuleExpression":
        coefficient = _fraction(coefficient)
        return RuleExpression(
            {key: coefficient * value for key, value in self.coefficients.items()},
            self.representatives,
        )

    def primitive(self) -> tuple["RuleExpression", Fraction]:
        """Return ``(primitive, scale)`` with ``self = scale * primitive``."""

        if not self.coefficients:
            raise ValueError("the zero expression has no primitive direction")
        from mfobs.core.linear_closure import primitive
        vector, scale = primitive(self.coefficients)
        return RuleExpression(vector, self.representatives), scale

    def record(self) -> dict:
        return {
            "term_count": len(self.coefficients),
            "terms": [
                {
                    "term_id": _term_id(key),
                    "coefficient": _encoded_coefficient(value),
                }
                for key, value in sorted(self.coefficients.items())
            ],
        }


def bilinear_commutator(
    observable: RuleExpression,
    generator: RuleExpression,
    pair_commutator: Callable[[object, object], object],
) -> tuple[RuleExpression, list[dict]]:
    """Expand a bilinear rule commutator without merging rate channels."""

    coefficients: Counter = Counter()
    representatives = {}
    pair_records = []
    for observable_key, observable_coefficient in sorted(observable.coefficients.items()):
        for generator_key, generator_coefficient in sorted(generator.coefficients.items()):
            result = pair_commutator(
                observable.representatives[observable_key],
                generator.representatives[generator_key],
            )
            scale = observable_coefficient * generator_coefficient
            for key, value in result.coefficients.items():
                coefficients[key] += scale * value
                representatives.setdefault(key, result.representatives[key])
            pair_records.append(
                {
                    "observable_term_id": _term_id(observable_key),
                    "generator_term_id": _term_id(generator_key),
                    "input_coefficient": _encoded_coefficient(scale),
                    "left": result.left.record(),
                    "right": result.right.record(),
                }
            )
    return RuleExpression(coefficients, representatives), pair_records


def jump_close_expression(
    expression: RuleExpression,
    *,
    jump_closure: Callable[[object], object] = sqpo_jump_closure,
) -> tuple[RuleExpression, list[dict]]:
    """Apply the selected jump closure termwise, retaining provenance."""

    coefficients: Counter = Counter()
    representatives = {}
    trace = []
    for source_key, coefficient in sorted(expression.coefficients.items()):
        jump = jump_closure(expression.representatives[source_key])
        target_key = jump.canonKey
        coefficients[target_key] += coefficient
        representatives.setdefault(target_key, jump)
        trace.append(
            {
                "raw_term_id": _term_id(source_key),
                "jump_term_id": _term_id(target_key),
                "coefficient": _encoded_coefficient(coefficient),
            }
        )
    return RuleExpression(coefficients, representatives), trace


@dataclass(frozen=True)
class ObservableDesignation:
    name: str
    meaning: str = ""
    metadata: Mapping | None = None


@dataclass
class ObservableClosureResult:
    status: str
    basis: dict[str, RuleExpression]
    origins: dict[str, dict]
    steps: list[dict]
    pending: list[str]
    evidence: list[dict]
    jump_projection: str
    obligations: list[dict] = field(default_factory=list)

    def record(self) -> dict:
        return {
            "status": self.status,
            "algorithm": "channel-residual-linear-closure-v1",
            "jump_projection": self.jump_projection,
            "arithmetic": "exact rational sparse rule-expression span",
            "basis": [
                {
                    "name": name,
                    "expression": expression.record(),
                    **self.origins[name],
                }
                for name, expression in self.basis.items()
            ],
            "steps": self.steps,
            "pending": self.pending,
            "obligations": self.obligations,
        }


def discover_observables(
    *,
    seeds: Mapping[str, RuleExpression],
    generator: RuleExpression,
    pair_commutator: Callable[[object, object], object],
    classify: Callable[[RuleExpression, int], ObservableDesignation] | None = None,
    jump_closure: Callable[[object], object] = sqpo_jump_closure,
    jump_projection_name: str | None = None,
    maximum_dimension: int = 64,
) -> ObservableClosureResult:
    """Single-channel specialization of Algorithm 1."""
    return discover_observables_by_channel(
        seeds=seeds, generators={"G": generator}, pair_commutator=pair_commutator,
        classify=classify, jump_closure=jump_closure,
        jump_projection_name=jump_projection_name, maximum_dimension=maximum_dimension,
    )



def discover_observables_by_channel(
    *,
    seeds: Mapping[str, RuleExpression],
    generators: Mapping[str, RuleExpression],
    pair_commutator: Callable[[object, object], object],
    classify: Callable[[RuleExpression, int], ObservableDesignation] | None = None,
    jump_closure: Callable[[object], object] = sqpo_jump_closure,
    jump_projection_name: str | None = None,
    maximum_dimension: int = 64,
) -> ObservableClosureResult:
    """Graph-expression backend for the shared channel-residual algorithm."""
    from mfobs.core.linear_closure import close_linear_span

    projection = jump_projection_name or getattr(jump_closure, "__name__", type(jump_closure).__name__)
    expressions = dict(seeds)
    representatives = {k: r for expr in seeds.values() for k, r in expr.representatives.items()}
    designations = {}
    evidence = []

    def image(source, channel):
        raw, pairs = bilinear_commutator(expressions[source], generators[channel], pair_commutator)
        jumped, trace = jump_close_expression(raw, jump_closure=jump_closure)
        representatives.update(jumped.representatives)
        evidence.append({"source": source, "channel": channel, "raw": raw, "jump_closed": jumped})
        return jumped.coefficients, {
            "operation": "p([source,G_channel])", "raw_commutator": raw.record(),
            "jump_closed_expression": jumped.record(), "composition_pairs": pairs,
            "jump_closure_trace": trace,
        }

    def expression(vector):
        return RuleExpression(vector, {k: representatives[k] for k in vector})

    def designate(vector, index):
        designation = classify(expression(vector), index) if classify else ObservableDesignation(f"O{index}")
        designations[designation.name] = designation
        return designation.name

    def added(name, vector, source, channel, before, scale):
        expressions[name] = expression(vector)

    result = close_linear_span(
        seeds={n: e.coefficients for n, e in seeds.items()}, channels=generators,
        image=image, encode=lambda v: expression(v).record(), designate=designate,
        on_add=added, maximum_dimension=maximum_dimension,
    )
    for name, origin in result.origins.items():
        designation = designations.get(name, ObservableDesignation(name))
        origin.update(meaning=designation.meaning, metadata=dict(designation.metadata or {}))
    return ObservableClosureResult(result.status, expressions, result.origins,
                                   result.steps, result.pending, evidence, projection,
                                   result.obligations)
