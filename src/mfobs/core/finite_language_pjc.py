"""Exact polynomial jump closure on a declared finite state-language quotient.

This backend is model-agnostic.  It receives exact observable functions and
rate-separated backward generators, discovers their unital pointwise
function algebra, groups transitions by observable jump, and interpolates each
jump propensity in that algebra.  No event delta or published equation is an
input.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Hashable, Mapping, Sequence

from mfobs.core.exact_span import ExactIncrementalSpan
from mfobs.core.polynomial_jump_closure import (
    JumpEvent,
    JumpModelProvenance,
    PolynomialJumpModel,
    PolynomialRing,
)


Vector = tuple[Fraction, ...]
Matrix = tuple[tuple[Fraction, ...], ...]
Monomial = tuple[int, ...]


def _pointwise(left: Vector, right: Vector) -> Vector:
    return tuple(a * b for a, b in zip(left, right))


def _pivot_columns(rows: Sequence[Sequence[Fraction]]) -> tuple[int, ...]:
    work = [list(row) for row in rows]
    if not work:
        return ()
    row_count, column_count = len(work), len(work[0])
    pivot_row = 0
    pivots = []
    for column in range(column_count):
        pivot = next((row for row in range(pivot_row, row_count)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return tuple(pivots)


def _inverse(rows: Matrix) -> Matrix:
    size = len(rows)
    if any(len(row) != size for row in rows):
        raise ValueError("inverse requires a square matrix")
    work = [list(row) + [Fraction(int(i == j)) for j in range(size)]
            for i, row in enumerate(rows)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


class _LinearBasis:
    def __init__(self, columns: Sequence[Vector]):
        self.columns = tuple(columns)
        if not self.columns:
            raise ValueError("basis must be nonempty")
        state_count = len(self.columns[0])
        if any(len(column) != state_count for column in self.columns):
            raise ValueError("basis columns have different lengths")
        transpose = tuple(tuple(column[row] for row in range(state_count))
                          for column in self.columns)
        self.pivot_rows = _pivot_columns(transpose)
        if len(self.pivot_rows) != len(self.columns):
            raise ValueError("basis columns are linearly dependent")
        square = tuple(tuple(column[row] for column in self.columns)
                       for row in self.pivot_rows)
        self.inverse = _inverse(square)

    def solve(self, target: Vector) -> tuple[Fraction, ...] | None:
        if len(target) != len(self.columns[0]):
            raise ValueError("target has the wrong state dimension")
        selected = tuple(target[row] for row in self.pivot_rows)
        coefficients = tuple(
            sum((entry * value for entry, value in zip(row, selected)), Fraction())
            for row in self.inverse
        )
        reconstructed = tuple(
            sum((coefficients[j] * self.columns[j][i]
                 for j in range(len(self.columns))), Fraction())
            for i in range(len(target))
        )
        return coefficients if reconstructed == target else None


@dataclass(frozen=True)
class FiniteFunctionAlgebra:
    labels: tuple[str, ...]
    monomials: tuple[Monomial, ...]
    vectors: tuple[Vector, ...]


@dataclass(frozen=True)
class FiniteLanguageCompilation:
    model: PolynomialJumpModel
    evidence: Mapping[str, Any]


def _monomial_label(monomial: Monomial, names: Sequence[str]) -> str:
    pieces = []
    for name, exponent in zip(names, monomial):
        if exponent == 1:
            pieces.append(name)
        elif exponent:
            pieces.append(f"{name}^{exponent}")
    return "*".join(pieces) or "1"


def discover_finite_function_algebra(
    names: tuple[str, ...], vectors: tuple[Vector, ...]
) -> FiniteFunctionAlgebra:
    """Discover a monomial basis of the generated unital function algebra."""

    if not names or len(names) != len(vectors) or len(names) != len(set(names)):
        raise ValueError("observable names and vectors must be nonempty and unique")
    span = ExactIncrementalSpan(vectors[0])
    for vector in vectors[1:]:
        if span.solve(vector) is not None:
            raise ValueError("basis columns are linearly dependent")
        span.add(vector)
    count = len(names)
    monomials = [tuple(int(i == j) for i in range(count)) for j in range(count)]
    columns = list(vectors)
    ones = tuple(Fraction(1) for _ in vectors[0])
    if span.solve(ones) is None:
        monomials.append((0,) * count)
        columns.append(ones)
        span.add(ones)
    cursor = 0
    while cursor < len(columns):
        for generator_index, generator in enumerate(vectors):
            candidate = _pointwise(columns[cursor], generator)
            if span.solve(candidate) is not None:
                continue
            exponent = list(monomials[cursor])
            exponent[generator_index] += 1
            monomials.append(tuple(exponent))
            columns.append(candidate)
            span.add(candidate)
            if len(columns) > len(ones):
                raise AssertionError("function algebra exceeded state-space dimension")
        cursor += 1
    return FiniteFunctionAlgebra(
        tuple(_monomial_label(value, names) for value in monomials),
        tuple(monomials),
        tuple(columns),
    )


def _signature_groups(vectors: tuple[Vector, ...]) -> tuple[tuple[int, ...], ...]:
    grouped = {}
    for state in range(len(vectors[0])):
        signature = tuple(vector[state] for vector in vectors)
        grouped.setdefault(signature, []).append(state)
    return tuple(tuple(group) for group in grouped.values())


def _lumpability(generator: Matrix, groups: tuple[tuple[int, ...], ...]) -> bool:
    for sources in groups:
        reference = tuple(sum(generator[sources[0]][target] for target in targets)
                          for targets in groups)
        for source in sources[1:]:
            actual = tuple(sum(generator[source][target] for target in targets)
                           for targets in groups)
            if actual != reference:
                return False
    return True


def compile_finite_language_jump_model(
    *,
    name: str,
    states: tuple[Hashable, ...],
    observable_names: tuple[str, ...],
    observable_vectors: tuple[Vector, ...],
    generators: Mapping[str, Matrix],
    input_definition_sha256: str = "",
    quotient_certificate: Mapping[str, Any],
    source_rules: Mapping[str, str],
) -> FiniteLanguageCompilation:
    """Derive event jumps and polynomial propensities from exact generators."""

    state_count = len(states)
    if not state_count or any(len(vector) != state_count for vector in observable_vectors):
        raise ValueError("observable functions do not cover the finite language")
    for channel, generator in generators.items():
        if len(generator) != state_count or any(len(row) != state_count for row in generator):
            raise ValueError(f"{channel}: generator has the wrong shape")
        if any(sum(row, Fraction()) for row in generator):
            raise ValueError(f"{channel}: generator rows must sum to zero")
        if any(generator[source][target] < 0
               for source in range(state_count) for target in range(state_count)
               if source != target):
            raise ValueError(f"{channel}: generator has a negative off-diagonal rate")
    if set(source_rules) != set(generators):
        raise ValueError("every generator channel needs one source rule")

    algebra = discover_finite_function_algebra(observable_names, observable_vectors)
    solver = _LinearBasis(algebra.vectors)
    groups = _signature_groups(observable_vectors)
    lumpability = {channel: _lumpability(generator, groups)
                   for channel, generator in generators.items()}
    if not all(lumpability.values()):
        raise ValueError("observable coordinate quotient is not strongly lumpable")

    ring = PolynomialRing(observable_names)
    events = []
    event_records = []
    for channel, generator in generators.items():
        propensities: dict[Vector, list[Fraction]] = {}
        for source in range(state_count):
            for target, rate in enumerate(generator[source]):
                if source == target or not rate:
                    continue
                delta = tuple(vector[target] - vector[source]
                              for vector in observable_vectors)
                if not any(delta):
                    continue
                propensities.setdefault(delta, [Fraction() for _ in states])[source] += rate
        for number, delta in enumerate(sorted(propensities), start=1):
            if any(value.denominator != 1 for value in delta):
                raise ValueError("pJC coordinates must have integral event jumps")
            propensity_vector = tuple(propensities[delta])
            coefficients = solver.solve(propensity_vector)
            if coefficients is None:
                raise ValueError("jump propensity is outside the generated function algebra")
            terms = {
                monomial: coefficient
                for monomial, coefficient in zip(algebra.monomials, coefficients)
                if coefficient
            }
            propensity = ring.polynomial(terms)
            identifier = f"{channel}-jump-{number}"
            events.append(JumpEvent(
                identifier, channel, tuple(int(value) for value in delta),
                propensity, source_rules[channel],
                "generic grouping of exact rule applications by coordinate jump",
            ))
            event_records.append({
                "event_class": identifier,
                "channel": channel,
                "delta": {observable: str(value)
                          for observable, value in zip(observable_names, delta) if value},
                "propensity": ring.record(propensity),
            })

    certificate = {
        **dict(quotient_certificate),
        "status": "VERIFIED_ALL_STATES",
        "finite_quotient_states": state_count,
        "observable_signature_classes": len(groups),
        "strong_lumpability_by_channel": lumpability,
        "function_algebra_dimension": len(algebra.vectors),
        "event_partition": "exact grouping by rule-application coordinate delta",
        "supplied_delta_table": False,
    }
    model = PolynomialJumpModel(
        name,
        ring,
        tuple(events),
        certificate,
        input_definition_sha256,
        JumpModelProvenance(
            backend="generic exact finite-language function algebra",
            event_partition_origin="exact grouping of literal rule applications by coordinate jump",
            jump_origin="target-minus-source evaluation after literal conditional-DPO application",
            activity_origin="exact finite function-algebra interpolation of rule-match multiplicities",
            jump_projection="conditional-DPO application on the declared finite quotient",
            composition_scope="linear discovery uses composition; pJC jumps use literal rule application",
            all_jumps_composition_derived=False,
            supplied_delta_table=False,
        ),
    )
    return FiniteLanguageCompilation(
        model,
        {
            "schema_version": 1,
            "status": "PASS",
            "backend": "generic exact finite-language function algebra",
            "input_definition_sha256": input_definition_sha256,
            "state_count": state_count,
            "observable_coordinates": list(observable_names),
            "observable_signature_classes": len(groups),
            "strong_lumpability_by_channel": lumpability,
            "generated_algebra": {
                "dimension": len(algebra.vectors),
                "basis": [
                    {"label": label, "degree": sum(monomial),
                     "monomial": {observable: exponent
                                  for observable, exponent in zip(observable_names, monomial)
                                  if exponent}}
                    for label, monomial in zip(algebra.labels, algebra.monomials)
                ],
            },
            "event_classes": event_records,
            "supplied_transition_or_delta_table": False,
        },
    )
