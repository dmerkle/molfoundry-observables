"""Cross-example checks for the shared linear/pJC architecture."""

import mfobs.models.remy.analysis as remy_analysis

from mfobs.core.polynomial_jump_closure import generator_action
from mfobs.models import brusselator_moments
from mfobs.models import brusselator_pjc
from mfobs.models import eight_slots_common
from mfobs.models import ring8_common
from mfobs.models import structural_nac_common
from mfobs.models.brusselator_growth_projection import BrusselatorGrowthGraphModel


def test_remy_uses_same_contract_for_both_branches():
    result = remy_analysis.common_analysis()
    digest = result["contract"]["definition_sha256"]
    assert result["linear"]["input_definition_sha256"] == digest
    assert result["pjc"]["input_definition_sha256"] == digest
    assert result["linear"]["status"] == "CLOSED"
    assert result["linear"]["basis"] == ["P3", "P2", "P1", "E"]
    assert result["pjc"]["operator"]["maximum_differential_order"] == 1
    assert result["pjc"]["operator_span"]["dimension"] == 4
    assert result["composition"]["jump_projection"] == "DPO with valid-tree input-count reduction"
    assert not result["composition"]["full_event_partition_composition_derived"]
    provenance = result["pjc"]["jump_model_provenance"]
    assert "E/P1" in provenance["composition_scope"]
    assert provenance["supplied_delta_table"] is False


def test_brusselator_pjc_closes_while_linear_moments_remain_open():
    result = brusselator_pjc.analyze(linear_limit=24)
    assert result["contract"]["analysis_branches"] == {
        "linear_observable_discovery": {
            "seeds": ["X", "Y"],
            "operation": "first commutator followed by diagonal jump projection",
        },
        "polynomial_jump_closure_emgf": {
            "targets": ["X", "Y"],
            "operation": "all nested diagonal commutators via polynomial event propensities",
        },
    }
    assert {
        observable["name"]: observable["role"]
        for observable in result["contract"]["observables"]
    } == {"X": "selected", "Y": "selected", "M": "support", "N": "support"}
    assert result["linear"]["status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert not {"M", "N"} & set(result["linear"]["basis"])
    assert result["pjc"]["status"] == "CLOSED"
    assert result["pjc"]["observable_coordinates"] == ["X", "Y", "M", "N"]
    assert result["pjc"]["event_class_count"] == 6
    assert result["pjc"]["operator"]["maximum_differential_order"] == 3
    assert result["composition"]["jump_projection"] == "DPO"
    assert result["pjc"]["jump_model_provenance"] == {
        "schema_version": 1,
        "backend": "MolFoundry DPO common-overlap composition",
        "event_partition_origin": (
            "application-supplied all-state structural match partition"
        ),
        "jump_origin": (
            "exact activity-identity coefficient of each composed observable commutator"
        ),
        "activity_origin": "application-supplied exact state-language reducer",
        "jump_projection": "DPO",
        "composition_scope": "all event-coordinate jumps",
        "all_jumps_composition_derived": True,
        "supplied_delta_table": False,
    }
    from mfobs.checks.brusselator import semantic_oracle_model
    from mfobs.core.composition_jump_model import compare_polynomial_jump_models
    assert compare_polynomial_jump_models(
        brusselator_pjc.polynomial_jump_model(), semantic_oracle_model())["status"] == "PASS"
    assert "independent_population_oracle" not in result["composition"]


def test_brusselator_structural_projection_reads_literal_component_inventories():
    graph_model = BrusselatorGrowthGraphModel(
        input_definition_sha256=brusselator_pjc.analysis_contract().digest
    )
    channels = {
        value.identifier: (dict(value.left), dict(value.right))
        for value in graph_model.population_channels()
    }
    assert channels == {
        "source": ({"A": 1}, {"X": 1}),
        "autocatalysis": ({"X": 2, "Y": 1}, {"X": 3}),
        "conversion": ({"B": 1, "X": 1}, {"D": 1, "Y": 1}),
        "sink": ({"X": 1}, {"E": 1}),
    }
    projection = graph_model.projection_record()
    assert projection["supplied_delta_table"] is False
    assert projection["population"]["construction"].endswith(
        "no observable jump vector is an input."
    )


def test_brusselator_first_generator_rows_match_existing_exact_adapter():
    model = brusselator_pjc.polynomial_jump_model()
    ring = model.ring
    for coordinate, moment in (("X", (1, 0)), ("Y", (0, 1))):
        for old_channel, new_channel in (
            ("source", "k1*a*V"),
            ("autocatalysis", "k2/V^2"),
            ("conversion", "k3*b"),
            ("sink", "k4"),
        ):
            actual = generator_action(model, ring.coordinate(coordinate), channel=new_channel)
            # Compare on a grid; the older adapter expresses the same exact
            # polynomial in the falling-factorial basis.
            expected_form = brusselator_moments.channel_derivative(moment, old_channel)
            for x in range(4):
                for y in range(4):
                    expected = sum(
                        coefficient
                        * brusselator_moments.factorial_expression(target).subs(
                            {brusselator_moments.NX: x, brusselator_moments.NY: y}
                        )
                        for target, coefficient in expected_form.items()
                    )
                    # Rate-channel parameters have been separated; the
                    # autocatalysis polynomial retains no V^-2 here.
                    if old_channel == "source":
                        expected = expected / (
                            brusselator_moments.A_FEED * brusselator_moments.VOLUME
                        )
                    elif old_channel == "autocatalysis":
                        expected = expected * brusselator_moments.VOLUME**2
                    elif old_channel == "conversion":
                        expected = expected / brusselator_moments.B_FEED
                    assert ring.evaluate(actual, {"X": x, "Y": y, "M": 0, "N": 0}) == expected


def test_ring8_uses_production_finite_rule_and_pjc_backends(ring_composition):
    result = ring8_common.analyze(composition_result=ring_composition)
    digest = result["contract"]["definition_sha256"]
    assert result["contract"]["analysis_branches"] == {
        "linear_observable_discovery": {
            "seeds": ["P0"],
            "operation": "first commutator followed by diagonal jump projection",
        },
        "polynomial_jump_closure_emgf": {
            "targets": ["P0"],
            "operation": "all nested diagonal commutators via polynomial event propensities",
        },
    }
    assert result["linear"]["input_definition_sha256"] == digest
    assert result["pjc"]["input_definition_sha256"] == digest
    assert result["linear"]["basis"][0] == "P0"
    assert result["linear"]["dimension"] == 30
    assert result["pjc"]["event_class_count"] == 196
    assert result["pjc"]["operator"]["maximum_differential_order"] == 1
    assert result["pjc"]["operator_span"]["dimension"] == 144
    assert (result["pjc"]["finite_language_compilation"]
            ["observable_signature_classes"] == 30)
    assert result["finite_rule_projection"]["supplied_transition_or_delta_table"] is False
    assert result["pjc"]["jump_model_provenance"]["supplied_delta_table"] is False
    assert result["cross_branch_validation"] == {
        "status": "PASS",
        "exact_scalar_checks": 2700,
        "comparison": (
            "literal conditional-DPO finite projection equals the independently "
            "composition-derived linear generator on all 30 quotient states"
        ),
    }


def test_eight_slots_uses_the_same_generated_finite_projection_backend(eight_composition):
    result = eight_slots_common.analyze(composition_result=eight_composition)
    assert result["linear"]["status"] == "CLOSED"
    assert result["linear"]["dimension"] == 48
    assert result["linear"]["semantic_rank"] == 47
    assert result["pjc"]["status"] == "CLOSED"
    assert result["pjc"]["event_class_count"] == 247
    assert result["pjc"]["finite_language_compilation"][
        "observable_signature_classes"
    ] == 83
    assert result["pjc"]["operator"]["maximum_differential_order"] == 1
    assert result["finite_rule_projection"][
        "supplied_transition_or_delta_table"
    ] is False
    assert result["cross_branch_validation"]["status"] == "PASS"


def test_growing_ring_uses_the_contract_and_reports_uncertified_closure():
    from mfobs import GrowingRing, analyze
    backend = GrowingRing()
    result = analyze(backend.specification(), backend=backend, worklist_limit=12).record()
    digest = result["contract"]["definition_sha256"]
    assert result["linear"]["input_definition_sha256"] == digest
    assert result["pjc"]["input_definition_sha256"] == digest
    assert result["linear"]["status"] == "UNKNOWN_EXPANSION_LIMIT"
    assert result["pjc"]["status"] == "UNKNOWN_UNBOUNDED_COORDINATE_LIFT"
    assert result["pjc"]["attempt_provenance"]["supplied_delta_table"] is False
    assert result["pjc"]["selected_coordinate_non_lumpability"]["status"] == "WITNESS"


def test_eight_slot_vectors_evaluate_signed_count_expressions(eight_composition):
    states = eight_slots_common.quotient_states()[:5]
    gml = eight_composition["observable_expressions"]["P0"][0]["gml"]
    one = [{"coefficient": "1", "gml": gml}]
    combined = [{"coefficient": "3/2", "gml": gml},
                {"coefficient": "-1/2", "gml": gml}]
    values = eight_slots_common._observable_vectors(
        {"observable_expressions": {"Q": one, "R": combined}}, states)
    assert values["Q"]==values["R"]
