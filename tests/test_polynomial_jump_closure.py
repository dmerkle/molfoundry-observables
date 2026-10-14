"""Model-independent checks for the production pJC/EMGF implementation."""
from fractions import Fraction

from mfobs.core.polynomial_jump_closure import (
    JumpEvent,
    JumpModelProvenance,
    PolynomialJumpModel,
    PolynomialRing,
    one_step_event_expansion,
)


def test_sparse_polynomial_ring_reduces_boolean_suffix_identities():
    ring = PolynomialRing(
        ("M", "S_A", "S_AA", "S_AB", "X"),
        boolean_variables=("S_A", "S_AA", "S_AB", "X"),
        exclusive_groups=(("S_AA", "S_AB"),),
        boolean_products={
            ("S_A", "S_AA"): "S_AA",
            ("S_A", "S_AB"): None,
        },
    )
    assert ring.text(ring.multiply(ring.coordinate("S_A"), ring.coordinate("S_A"))) == "S_A"
    assert ring.text(ring.multiply(ring.coordinate("S_A"), ring.coordinate("S_AA"))) == "S_AA"
    assert ring.multiply(ring.coordinate("S_A"), ring.coordinate("S_AB")).is_zero
    assert ring.text(ring.multiply(
        ring.coordinate("M"), ring.coordinate("S_AA"), ring.coordinate("X")
    )) == "M*S_AA*X"


def test_polynomial_shift_is_exact_before_boolean_reduction():
    ring = PolynomialRing(("M", "X"), boolean_variables=("X",))
    value = ring.multiply(ring.coordinate("M"), ring.coordinate("X"))
    shifted = ring.shift(value, (-1, 0))
    assert ring.text(shifted) == "-X + M*X"
    assert ring.evaluate(shifted, {"M": 7, "X": 1}) == Fraction(6)


def test_generic_one_step_expansion_evaluates_activities_and_jumps():
    ring = PolynomialRing(("X",))
    model = PolynomialJumpModel(
        "birth-death",
        ring,
        (
            JumpEvent("birth", "alpha", (1,), ring.constant(1), "birth", "always"),
            JumpEvent("death", "beta", (-1,), ring.coordinate("X"), "death", "one X"),
        ),
        {"status": "VERIFIED_ALL_STATES"},
        "test-input",
        JumpModelProvenance(
            "test", "test", "test", "test", "DPO", "test", True, False
        ),
    )
    result = one_step_event_expansion(model, {"X": 2})
    assert [(row["event_class"], row["activity"], row["target"]) for row in result["transitions"]] == [
        ("birth", "1", {"X": 3}),
        ("death", "2", {"X": 1}),
    ]
