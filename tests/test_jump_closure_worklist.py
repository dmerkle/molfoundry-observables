"""Model-independent exact-span checks for first-commutator discovery."""
from fractions import Fraction
from mfobs.core.linear_closure import RationalSpan

from mfobs.core.jump_closure_worklist import (
    RuleExpression,
)


def expression(coefficients):
    return RuleExpression(coefficients, {key: object() for key in coefficients})


def test_primitive_rule_expression_uses_exact_rational_arithmetic():
    value = expression({"a": Fraction(-2, 3), "b": Fraction(-4, 3)})
    primitive, scale = value.primitive()
    assert primitive.coefficients == {"a": Fraction(1), "b": Fraction(2)}
    assert scale == Fraction(-2, 3)


def test_sparse_span_reduces_exactly_when_a_later_pivot_sorts_first():
    span = RationalSpan()
    first = expression({"m": 1, "z": 1})
    second = expression({"a": 1, "z": 1})
    assert span.insert("first", first.coefficients)
    assert span.insert("second", second.coefficients)
    target = expression({"a": 3, "m": -2, "z": 1})
    residual, coordinates = span.reduce(target.coefficients)
    assert not residual
    assert coordinates == {"first": Fraction(-2), "second": Fraction(3)}
