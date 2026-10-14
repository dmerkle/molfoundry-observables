"""Generic finite-rule and finite-function pJC backend checks."""
from mfobs.core.graph import node_rows, make_graph
from mfobs.core.rules import make_rule
from fractions import Fraction


from mfobs.core.finite_language_pjc import compile_finite_language_jump_model
from mfobs.core.finite_rule_projection import (
    FiniteRuleChannel,
    project_finite_rule_generators,
)

from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure


def test_literal_rule_application_derives_generator_and_pjc_without_delta_input():
    token = make_graph({0: "X"})
    remove = make_rule("remove", token, make_graph({}))

    projection = project_finite_rule_generators(
        states=(0, 1),
        channels=(FiniteRuleChannel("remove", remove, Fraction(1)),),
        representative=lambda state, _channel: make_graph(
            {index: "X" for index in range(state)}
        ),
        project_target=lambda _source, _channel, target: len(node_rows(target)),
    )
    assert projection.generators["remove"] == (
        (Fraction(0), Fraction(0)),
        (Fraction(1), Fraction(-1)),
    )
    compilation = compile_finite_language_jump_model(
        name="two-state-remove",
        states=(0, 1),
        observable_names=("X",),
        observable_vectors=((Fraction(0), Fraction(1)),),
        generators=projection.generators,
        input_definition_sha256="test",
        quotient_certificate={"method": "complete two-state language"},
        source_rules={"remove": "remove"},
    )
    report = analyze_polynomial_jump_closure(compilation.model)
    assert report["status"] == "CLOSED"
    assert report["event_class_count"] == 1
    assert report["operator"]["terms"][0]["delta"] == {"X": -1}
    assert report["operator"]["maximum_differential_order"] == 1
    assert not compilation.evidence["supplied_transition_or_delta_table"]


def test_finite_quotient_with_zero_generators_has_a_closed_empty_event_model():
    compilation = compile_finite_language_jump_model(
        name="stationary", states=(0, 1), observable_names=("X",),
        observable_vectors=((Fraction(0), Fraction(1)),),
        generators={"k": ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))},
        quotient_certificate={"method": "complete two-state language"},
        source_rules={"k": "preserve"},
    )
    report = analyze_polynomial_jump_closure(compilation.model)
    assert report["status"] == "CLOSED"
    assert report["event_class_count"] == report["operator_span"]["dimension"] == 0
    assert report["first_moment_generator_images"] == {"X": {}}
