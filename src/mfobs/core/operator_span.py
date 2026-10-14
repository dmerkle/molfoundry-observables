"""Exact worklists for commutator-closed finite operator spans.

The algorithm is independent of graphs and of any particular model. A model
supplies a finite semantic event partition and the jump of each requested
diagonal observable on every event class. On this quotient, ``ad_O`` acts by
pointwise multiplication with the jump vector. The FIFO worklist discovers
the generated operator span and exact adjoint matrices; it does not receive a
target basis or published commutator identities.

For graph-rule models the quotient needs a separate, all-state certificate.
Finite observations are useful regression tests, but are not such a
certificate by themselves.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Mapping, Sequence

from mfobs.core.exact_span import ExactIncrementalSpan


Vector = tuple[Fraction, ...]


def _vector(values: Iterable[int | Fraction]) -> Vector:
    return tuple(Fraction(value) for value in values)


def _encode_vector(values: Sequence[Fraction]) -> list[str]:
    return [str(value) for value in values]


def _solve_columns(columns: Sequence[Vector], target: Vector) -> tuple[Fraction, ...] | None:
    """Solve ``columns * x = target`` exactly, returning ``None`` if outside."""

    if not columns:
        return () if not any(target) else None
    rows = len(target)
    width = len(columns)
    matrix = [
        [columns[column][row] for column in range(width)] + [target[row]]
        for row in range(rows)
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, rows) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    if any(not any(row[:width]) and row[width] for row in matrix):
        return None
    solution = [Fraction(0) for _ in range(width)]
    for row, column in enumerate(pivot_columns):
        solution[column] = matrix[row][width]
    if any(
        sum(columns[column][row] * solution[column] for column in range(width))
        != target[row]
        for row in range(rows)
    ):
        return None
    return tuple(solution)


_IncrementalSpan = ExactIncrementalSpan


@dataclass(frozen=True)
class OperatorBasisElement:
    """One worklist-discovered operator before jump closure."""

    name: str
    vector: Vector
    word: tuple[str, ...]
    source: str | None
    observable: str | None

    @property
    def expression(self) -> str:
        value = "G"
        for observable in reversed(self.word):
            value = f"ad_{observable}({value})"
        return value


@dataclass(frozen=True)
class OperatorSpanResult:
    event_classes: tuple[str, ...]
    observables: tuple[str, ...]
    basis: tuple[OperatorBasisElement, ...]
    adjoint_matrices: Mapping[str, tuple[tuple[Fraction, ...], ...]]
    trace: tuple[dict, ...]

    def coordinates(self, vector: Sequence[int | Fraction]) -> tuple[Fraction, ...]:
        result = _solve_columns(
            [element.vector for element in self.basis], _vector(vector)
        )
        if result is None:
            raise ValueError("operator is outside the discovered span")
        return result

    def record(self) -> dict:
        return {
            "status": "CLOSED",
            "algorithm": "exact FIFO span worklist under diagonal commutator actions",
            "event_classes": list(self.event_classes),
            "observables": list(self.observables),
            "dimension": len(self.basis),
            "basis": [
                {
                    "name": element.name,
                    "expression": element.expression,
                    "commutator_word_outer_to_inner": list(element.word),
                    "event_coefficients": _encode_vector(element.vector),
                    "raw_rule_expression_before_jump_closure": [
                        {"event_rule": event_class, "coefficient": str(coefficient)}
                        for event_class, coefficient in zip(
                            self.event_classes, element.vector
                        )
                        if coefficient
                    ],
                    "source": element.source,
                    "observable": element.observable,
                }
                for element in self.basis
            ],
            "adjoint_matrices": {
                name: [[str(value) for value in row] for row in matrix]
                for name, matrix in self.adjoint_matrices.items()
            },
            "trace": list(self.trace),
        }


def discover_diagonal_commutator_span(
    *,
    event_classes: Sequence[str],
    jumps: Mapping[str, Sequence[int | Fraction]],
    seed: Sequence[int | Fraction] | None = None,
) -> OperatorSpanResult:
    """Discover the span generated from ``seed`` by every ``ad_O`` action.

    ``jumps[O][j]`` is the change of observable ``O`` in semantic event class
    ``j``. No closure dimension, basis vector, or commutator identity is an
    input to this routine.
    """

    classes = tuple(event_classes)
    observable_names = tuple(jumps)
    jump_vectors = {name: _vector(values) for name, values in jumps.items()}
    if any(len(values) != len(classes) for values in jump_vectors.values()):
        raise ValueError("every jump vector must cover every event class")
    initial = _vector(seed if seed is not None else [1] * len(classes))
    if not classes and not initial:
        return OperatorSpanResult(classes, observable_names, (),
                                  {name: () for name in observable_names}, ())
    if len(initial) != len(classes) or not any(initial):
        raise ValueError("the seed must be a nonzero event-class vector")

    basis = [OperatorBasisElement("B0", initial, (), None, None)]
    span = _IncrementalSpan(initial)
    pending = deque([0])
    trace: list[dict] = []
    while pending:
        source_index = pending.popleft()
        source = basis[source_index]
        for observable in observable_names:
            candidate = tuple(
                value * jump
                for value, jump in zip(source.vector, jump_vectors[observable])
            )
            coordinates = span.solve(candidate)
            record = {
                "source": source.name,
                "source_expression": source.expression,
                "observable": observable,
                "raw_event_coefficients": _encode_vector(candidate),
            }
            if coordinates is None:
                target = OperatorBasisElement(
                    f"B{len(basis)}",
                    candidate,
                    (observable, *source.word),
                    source.name,
                    observable,
                )
                basis.append(target)
                span.add(candidate)
                pending.append(len(basis) - 1)
                record.update({"decision": "ADD", "target": target.name})
            else:
                record.update(
                    {
                        "decision": "DEPENDENT",
                        "coordinates": _encode_vector(coordinates),
                    }
                )
            trace.append(record)

    matrices: dict[str, tuple[tuple[Fraction, ...], ...]] = {}
    for observable in observable_names:
        columns = []
        for element in basis:
            candidate = tuple(
                value * jump
                for value, jump in zip(element.vector, jump_vectors[observable])
            )
            coordinates = span.solve(candidate)
            if coordinates is None:
                raise AssertionError("worklist terminated before closure")
            columns.append(coordinates)
        matrices[observable] = tuple(
            tuple(columns[column][row] for column in range(len(basis)))
            for row in range(len(basis))
        )

    return OperatorSpanResult(
        classes, observable_names, tuple(basis), matrices, tuple(trace)
    )


def linear_propensity_record(
    coefficients: Sequence[int | Fraction], observable_names: Sequence[str]
) -> dict:
    values = _vector(coefficients)
    if len(values) != len(observable_names):
        raise ValueError("propensity coefficient dimension mismatch")
    terms = [
        {"observable": name, "coefficient": str(value)}
        for name, value in zip(observable_names, values)
        if value
    ]
    return {"terms": terms, "text": _linear_text(values, observable_names)}


def _linear_text(coefficients: Sequence[Fraction], names: Sequence[str]) -> str:
    parts: list[str] = []
    for coefficient, name in zip(coefficients, names):
        if not coefficient:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = abs(coefficient)
        body = name if magnitude == 1 else f"{magnitude}*{name}"
        if not parts:
            parts.append(body if sign == "+" else f"-{body}")
        else:
            parts.append(f" {sign} {body}")
    return "".join(parts) or "0"


def emgf_operator_record(
    *,
    observable_names: Sequence[str],
    event_classes: Sequence[Mapping],
    variables: Sequence[str],
    jump_projection: str = "SqPO jump closure",
) -> dict:
    """Apply jump closure and export weighted-growth and CTMC EMGF terms."""

    names = tuple(observable_names)
    variables = tuple(variables)
    if len(names) != len(variables):
        raise ValueError("one EMGF variable is required per observable")
    weighted = []
    ctmc = []
    for event in event_classes:
        delta = tuple(int(value) for value in event["delta"])
        coefficients = tuple(Fraction(value) for value in event["propensity"])
        if len(delta) != len(names) or len(coefficients) != len(names):
            raise ValueError("event-class dimension mismatch")
        exponent_terms = [
            (variable, value)
            for variable, value in zip(variables, delta)
            if value
        ]
        exponent = " + ".join(
            variable if value == 1 else f"-{variable}" if value == -1 else f"{value}*{variable}"
            for variable, value in exponent_terms
        ) or "0"
        factor = f"exp({exponent})"
        propensity = linear_propensity_record(coefficients, names)
        base = {
            "event_class": event["id"],
            "delta": list(delta),
            "exponential_factor": factor,
            "propensity": propensity,
        }
        weighted.append(base)
        ctmc.append({**base, "exponential_factor": f"({factor} - 1)"})

    derivative_coefficients = []
    for index, name in enumerate(names):
        contributions = []
        for term in weighted:
            coefficient = next(
                (
                    Fraction(item["coefficient"])
                    for item in term["propensity"]["terms"]
                    if item["observable"] == name
                ),
                Fraction(0),
            )
            if coefficient:
                contributions.append(
                    {
                        "event_class": term["event_class"],
                        "coefficient": str(coefficient),
                        "factor": term["exponential_factor"],
                    }
                )
        derivative_coefficients.append(
            {"observable": name, "variable": variables[index], "terms": contributions}
        )
    return {
        "construction": (
            f"{jump_projection} of the discovered event-class operator span; "
            "weighted growth keeps exp(delta·omega), CTMC uses exp(delta·omega)-1"
        ),
        "weighted_growth_terms": weighted,
        "ctmc_terms": ctmc,
        "expanded_derivative_coefficients": derivative_coefficients,
    }
