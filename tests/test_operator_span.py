"""Exact regression checks for the incremental operator-span basis."""
from fractions import Fraction
import random

from mfobs.core.operator_span import _IncrementalSpan, _solve_columns


def test_incremental_span_matches_reference_elimination():
    randomizer = random.Random(271828)
    for dimension in range(2, 9):
        columns = []
        span = None
        for _ in range(3 * dimension):
            candidate = tuple(
                Fraction(randomizer.randint(-3, 3)) for _ in range(dimension)
            )
            if not any(candidate):
                continue
            reference = _solve_columns(columns, candidate)
            if span is None:
                assert reference is None
                columns.append(candidate)
                span = _IncrementalSpan(candidate)
            elif reference is None:
                assert span.solve(candidate) is None
                columns.append(candidate)
                span.add(candidate)
            else:
                assert span.solve(candidate) == reference
            if len(columns) == dimension:
                break

        assert span is not None
        for _ in range(20):
            target = tuple(
                Fraction(randomizer.randint(-5, 5)) for _ in range(dimension)
            )
            assert span.solve(target) == _solve_columns(columns, target)
