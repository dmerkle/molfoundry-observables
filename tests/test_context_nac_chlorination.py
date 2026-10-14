"""Regression tests for the raw-composition context-NAC chlorination pipeline."""
from __future__ import annotations

import pytest

import mfobs
from mfobs.models.chlorination.analysis import Chlorination
from mfobs.checks.chlorination import validate_all_masks


@pytest.mark.parametrize("variant, guarded", [("guarded", True), ("ablated", False)])
def test_public_chlorination_selects_the_nac_model(variant, guarded):
    model = mfobs.Chlorination(variant)
    backend = Chlorination(variant="guarded" if guarded else "ablated")
    assert model.rule_gmls() == {
        "chlorinate": backend.rule.getGMLString(),
        "grow": backend.growth_rule.getGMLString(),
    }
    assert model.query_gmls() == {"Z": backend.patterns["Z"].rule.getGMLString()}
    assert "constrainNoEdge" not in model.rule_gmls()["chlorinate"]


def test_unknown_chlorination_variant_is_rejected():
    with pytest.raises(ValueError, match="unknown chlorination variant"):
        mfobs.Chlorination("typo")


@pytest.mark.parametrize("filters", [{}, {"nodeLabels": ["C", "C-"]}, {"edgeLabels": [":"]}])
def test_omitted_adjacency_filters_are_unsupported_instead_of_exact_zero(filters):
    import molfoundry as mf
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.core.pattern_support import UnsupportedReduction
    from mfobs.models.chlorination.inputs import molecule, seed_graph
    from mfobs.models.chlorination.counting import LocalConditionalReducer

    condition = mf.AdjacencyConstraint(0, mf.ConstraintOp(">="), 1, **filters)
    rule = mf.Rule.identity(seed_graph(), "ring_neighbour", constraints=[condition], add=False)
    assert rule.numMatches([molecule(0, 1)], mf.ApplyMode.UseAny) == 1
    reducer = LocalConditionalReducer()
    with pytest.raises(UnsupportedCondition, match="explicit.*filters"):
        reducer.reduce(rule)
    with pytest.raises(UnsupportedReduction, match="explicit.*filters"):
        reducer.reduce_image(rule)


@pytest.fixture(scope="module")
def backend():
    return mfobs.Chlorination()


@pytest.fixture(scope="module")
def analysis(backend):
    return mfobs.analyze(backend.specification(), backend=backend).record()


def test_context_nac_support_closes_from_actual_commutators(analysis):
    result = analysis

    assert result["status"] == "PASS"
    assert result["rule"]["condition_kind"] == "generic-not-exists"

    bridge = result["composition_bridge"]
    assert bridge["all_generator_images_from_raw_commutators"] is True
    assert bridge["all_projected_terms_reduced_from_engine_conditions"] is True
    assert bridge["activity_from_dpo_jump_closure"] is True
    assert bridge["state_space_enumeration"] is False

    support = result["synthesis"]["support"]
    assert support["status"] == "POLYNOMIAL_SUPPORT_CLOSED"
    assert list(support["patterns"]) == ["Z", "P1", "P2", "P3"]
    assert support["pending"] == []

    assert support["activities"]["chlorinate"]["text"] == (
        "-P1*Z + 4*P1 - P2*Z + 4*P2"
    )
    images = support["generator_images"]
    assert images["Z"]["chlorinate"]["text"] == "-P1*Z + 4*P1 - P2*Z + 4*P2"
    assert images["P1"]["chlorinate"]["text"] == "P1*Z - 4*P1 + P2*Z - 4*P2"
    assert images["P2"]["chlorinate"]["text"] == "2*P2*Z - 8*P2 + P3*Z - 4*P3"
    assert images["P3"]["chlorinate"]["text"] == "4*P3*Z - 16*P3"

    assert set(bridge["generator_reductions"]) == set(support["patterns"])
    assert all(row["projected_conditional_terms"] > 0
               for row in bridge["generator_reductions"].values())


def test_context_nac_events_recover_raw_composition_generator(analysis):
    result = analysis
    events = result["events"]

    assert {k: v for k, v in result["pjc"].items() if k != "joint_partition"} == {
        "status": "CLOSED_FINITE_EVENT_MODEL",
        "coordinates": ["Z", "P1", "P2", "P3"],
        "event_count": 4,
        "reason": None,
        "all_activities_polynomial": True,
        "constant_integer_jumps": True,
        "generator_match": "exact equality with raw-composition generator images",
    }
    proof = result["pjc"]["joint_partition"]
    assert proof["status"] == "VERIFIED_SYMBOLIC_JOINT_PARTITION"
    assert proof["disjoint"] and proof["exhaustive"]
    assert proof["state_space_enumeration"] is False
    assert {tuple(row["delta"].get(n, 0) for n in ("Z", "P1", "P2", "P3"))
            for row in events} == {
        (1, -1, 0, 0), (1, -1, -1, 0), (1, -1, -1, -2), (1, -1, -2, -2)}


def test_context_nac_post_synthesis_validation_is_independent_and_complete(backend, analysis):
    validation = validate_all_masks(backend, analysis)

    assert validation["status"] == "PASS"
    assert "post-synthesis validation only" in validation["role"]
    assert validation["masks_checked"] == 16
    assert validation["generator_scalar_checks"] == 192
    assert validation["molfoundry_reachability_checks"] == 48
    assert validation["growth_reachability_checks"] == 48
    assert validation["states_checked"] == 48
    assert validation["failures"] == []


@pytest.mark.parametrize("budgets", [dict(maximum_patterns=1), dict(maximum_patterns=3),
                                    dict(maximum_expansions=0), dict(maximum_expansions=3)])
def test_incomplete_support_has_no_closed_events(budgets):
    backend2 = Chlorination()
    result = mfobs.analyze(backend2.specification(), backend=backend2, **budgets).record()
    assert result["status"] == "UNKNOWN"
    assert result["synthesis"]["support"]["pending"]
    assert result["pjc"]["status"] == "UNKNOWN"
    assert result["pjc"]["joint_partition"] is None
    assert result["events"] == []
    assert "post_synthesis_validation" not in result


@pytest.mark.parametrize("successful_calls", [0, 1])
def test_unsupported_reduction_keeps_pending_work(monkeypatch, successful_calls):
    from mfobs.core.conditions import UnsupportedCondition
    backend = Chlorination()
    original = backend.reducer.reduce
    def unsupported(*args):
        nonlocal successful_calls
        if successful_calls:
            successful_calls -= 1
            return original(*args)
        raise UnsupportedCondition("outside the declared fragment")
    monkeypatch.setattr(backend.reducer, "reduce", unsupported)
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    assert result["synthesis"]["support"]["status"] == "UNKNOWN_UNSUPPORTED_REDUCTION"
    assert result["events"] == []


def test_synthesis_does_not_call_finite_validation_helpers(monkeypatch):
    def forbidden(*args):
        raise AssertionError("state enumeration was used by synthesis")
    from mfobs.checks import chlorination as check
    for name in ("_molecule", "_enabled_positions", "validate_all_masks"):
        monkeypatch.setattr(check, name, forbidden)
    from mfobs.checks import chlorination as checks
    monkeypatch.setattr(checks, "_coordinate_values", forbidden)
    from mfobs.models.chlorination import inputs as chemistry
    monkeypatch.setattr(chemistry, "molecule", forbidden)
    monkeypatch.setattr(chemistry, "host", forbidden)
    backend3 = Chlorination()
    result = mfobs.analyze(backend3.specification(), backend=backend3).record()
    assert result["pjc"]["joint_partition"]["status"] == "VERIFIED_SYMBOLIC_JOINT_PARTITION"


def test_equal_first_moments_do_not_certify_joint_events():
    from dataclasses import replace
    import sympy as sp
    from mfobs.models.chlorination.events import _derive_events, _certify_events
    backend = Chlorination()
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    names = result["pjc"]["coordinates"]
    good = _derive_events(backend.rule, backend.reducer, names)
    first = good[0]
    # Replace a jump by a double jump at half rate and a no-op at half rate.
    # Total activity and every generator image agree, but the law is wrong.
    bad = (replace(first, activity=first.activity / 2,
                   delta=tuple(2 * x for x in first.delta)),
           replace(first, activity=first.activity / 2, delta=(0, 0, 0, 0)), *good[1:])
    assert sp.expand(sum(e.activity for e in bad) - sum(e.activity for e in good)) == 0
    for i in range(4):
        assert sp.expand(sum(e.activity * e.delta[i] for e in bad)
                         - sum(e.activity * e.delta[i] for e in good)) == 0
    with pytest.raises(AssertionError, match="misses an enabled match"):
        _certify_events(backend.rule, bad, backend.reducer, names)


@pytest.mark.parametrize("relabel", [False, True])
def test_condition_indices_ignore_gml_ids_and_edge_order(relabel):
    import re
    import molfoundry as mf
    from mfobs.core.engine import dpo_jump_closure
    from mfobs.models.chlorination.inputs import RULE_GML
    from mfobs.models.chlorination.counting import _count_on_symbolic_host, RING
    text = RULE_GML.replace(
        'edge [ source 0 target 4 label ":" ]\n        edge [ source 4 target 5 label "-" ]',
        'edge [ source 4 target 5 label "-" ]\n        edge [ source 0 target 4 label ":" ]')
    if relabel:
        ids = [90, 12, 700, 3, 21, 44]
        text = re.sub(r'\b(id|source|target) (\d+)',
                      lambda m: m[1] + ' ' + str(ids[int(m[2])]), text)
    rule = mf.Rule.fromGMLString(text, add=False)
    value = _count_on_symbolic_host(dpo_jump_closure(rule), 1)
    assert value == RING.polynomial({(1,1,0,0): 1, (0,0,1,1): 1,
                                     (1,1,1,0): 1, (0,1,1,1): 1})
    from mfobs.models.chlorination.events import _derive_events, _certify_events
    backend = Chlorination()
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    names = result["pjc"]["coordinates"]
    events = _derive_events(rule, backend.reducer, names)
    assert _certify_events(rule, events, backend.reducer, names)["status"] == "VERIFIED_SYMBOLIC_JOINT_PARTITION"


def test_existential_wildcards_are_extension_relative():
    from molfoundry.constraints import ConditionEdge, ConditionExists, ConditionTrue, ConditionVertex
    from mfobs.models.chlorination.counting import _condition_value, ONE
    condition = ConditionExists(vertices=(ConditionVertex(id=99, label=None),),
                                edges=(ConditionEdge(source=0, target=99, label="-"),),
                                condition=ConditionTrue())
    assert _condition_value(condition, {0: 5}, {0}, {5: "C", 6: "Cl"}, {},
                            {(5, 6): "-"}) == ONE


def test_unanchored_existential_is_not_counted_on_a_truncated_host():
    import molfoundry as mf
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.models.chlorination.counting import LocalConditionalReducer
    rule = mf.Rule.fromGMLString('''rule [ context [ node [ id 0 label "C-" ] ]
        condition [ not [ exists [ node [ id 1 label "Cl" ] ] ] ] ]''', add=False)
    with pytest.raises(UnsupportedCondition, match="Unanchored existential"):
        LocalConditionalReducer().reduce(rule)


def test_nested_condition_extends_the_current_anchor():
    import molfoundry as mf
    from mfobs.core.engine import dpo_jump_closure
    from mfobs.models.chlorination.inputs import RULE_GML
    from mfobs.models.chlorination.counting import _count_on_symbolic_host, H, RING, ONE, ZERO
    text = RULE_GML.split("  condition [")[0] + '''condition [ not [ exists [
      node [ id 17 label "C" ] edge [ source 0 target 17 label ":" ]
      not [ exists [ node [ id 90 label "Cl" ]
                     edge [ source 17 target 90 label "-" ] ] ]
    ] ] ] ]'''
    rule = mf.Rule.fromGMLString(text, add=False)
    # Every neutral ring neighbour must now be chlorinated.
    expected = ZERO
    for i in range(4):
        site = H[i]
        for j in (i - 1, i + 1):
            if j in range(4):
                site = RING.multiply(site, RING.subtract(ONE, H[j]))
        expected = RING.add(expected, site)
    assert _count_on_symbolic_host(dpo_jump_closure(rule), 1) == expected


def test_analysis_does_not_execute_postchecks(monkeypatch):
    from mfobs.checks import chlorination as check
    def forbidden(*args):
        raise AssertionError("postcheck called during analysis")
    monkeypatch.setattr(check, "validate_all_masks", forbidden)
    backend4 = Chlorination()
    result = mfobs.analyze(backend4.specification(), backend=backend4).record()
    assert result["status"] == "PASS"
    assert "post_synthesis_validation" not in result


def test_only_the_query_is_supplied_and_positive_patterns_are_generated():
    from mfobs.models.chlorination.counting import RING
    backend = Chlorination()
    assert list(backend.patterns) == ["Z"]
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    assert result["status"] == "PASS"
    # Expected functions are an independent regression oracle, not model input.
    assert backend.reducer.coordinates["P1"] == RING.polynomial({(1,1,0,0): 1, (0,0,1,1): 1})
    assert backend.reducer.coordinates["P2"] == RING.polynomial({(1,1,1,0): 1, (0,1,1,1): 1})
    assert backend.reducer.coordinates["P3"] == RING.polynomial({(1,1,1,1): 2})
    for name in ("P1", "P2", "P3"):
        text = backend.patterns[name].rule.getGMLString()
        assert "condition [" not in text
        assert "constrain" not in text


def test_removing_nac_changes_support_and_events_without_changing_backend():
    backend = Chlorination("ablated")
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    assert result["status"] == "PASS"
    assert list(result["synthesis"]["support"]["patterns"]) == ["Z"]
    assert list(backend.patterns) == ["Z"]
    assert result["pjc"]["event_count"] == 1
    assert result["events"][0]["delta"] == {"Z": 1}
    assert result["events"][0]["activity"]["text"] == "Z**2 - 8*Z + 16"
    validation = validate_all_masks(backend, result)
    assert validation["status"] == "PASS"
    assert validation["generator_scalar_checks"] == 48


def test_canonical_pattern_generation_identifies_reflections():
    from mfobs.models.chlorination.counting import LocalConditionalReducer, RING
    reducer = LocalConditionalReducer()
    local = RING.polynomial({(1,1,0,0): 1, (0,0,1,1): 1})
    expr = reducer.express(local)
    assert len(reducer.patterns) == 2
    assert reducer.express(RING.scale(local, 3)) == 3 * expr
    assert len(reducer.patterns) == 2


def test_unrepresentable_joint_rates_do_not_claim_closure(monkeypatch):
    from mfobs.models.chlorination import analysis as pipeline
    from mfobs.core.conditions import UnsupportedCondition
    def unsupported(*args):
        raise UnsupportedCondition("Joint activities require additional positive patterns")
    monkeypatch.setattr(pipeline, "_derive_events", unsupported)
    backend5 = Chlorination()
    result = mfobs.analyze(backend5.specification(), backend=backend5).record()
    assert result["synthesis"]["support"]["status"] == "POLYNOMIAL_SUPPORT_CLOSED"
    assert result["status"] == result["pjc"]["status"] == "UNKNOWN"
    assert result["events"] == []
    assert "additional positive patterns" in result["pjc"]["reason"]
    assert "post_synthesis_validation" not in result


def test_event_transport_rejects_a_different_bond_exchange():
    import molfoundry as mf
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.models.chlorination.inputs import RULE_GML
    from mfobs.models.chlorination.events import _match_jumps
    backend = Chlorination()
    text = RULE_GML.replace('edge [ source 0 target 3 label "-" ]',
                            'edge [ source 0 target 2 label "-" ]')
    rule = mf.Rule.fromGMLString(text, add=False)
    with pytest.raises(UnsupportedCondition, match="Unsupported source bond exchange"):
        list(_match_jumps(rule, backend.reducer.coordinates))


def test_unbounded_language_and_active_growth_are_certified(analysis):
    language = analysis["synthesis"]["state_language"]
    assert language["finite"] is False
    assert "unbounded" in language["tail_extension"]
    growth = analysis["growth"]
    assert growth["certificate"]["status"] == "STRUCTURAL_ZERO_JUMP"
    assert growth["certificate"]["all_tail_lengths"]
    assert growth["generator_images"] == {n: "0" for n in analysis["pjc"]["coordinates"]}
    bridge = analysis["composition_bridge"]
    reductions = [bridge["activity_reduction"]] + [term["reduction"]
        for row in bridge["generator_reductions"].values() for term in row["terms"]]
    assert all(r["locality"]["status"] == "STRUCTURAL_ALL_TAIL_LENGTHS" for r in reductions)
    assert any(r["locality"]["checked_existentials"] > 0 for r in reductions)


def test_connected_nac_can_still_inspect_the_unbounded_tail():
    import molfoundry as mf
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.models.chlorination.counting import LocalConditionalReducer
    # The first tail carbon has three H at n=0, but only two at n>=1.
    # Connectedness to C- alone would incorrectly authorize stub evaluation.
    rule = mf.Rule.fromGMLString('''rule [ context [ node [ id 0 label "C-" ] ]
      condition [ not [ exists [
        node [ id 1 label "C" ] node [ id 2 label "H" ]
        node [ id 3 label "H" ] node [ id 4 label "H" ]
        edge [ source 0 target 1 label "-" ]
        edge [ source 1 target 2 label "-" ] edge [ source 1 target 3 label "-" ]
        edge [ source 1 target 4 label "-" ]
      ] ] ] ]''', add=False)
    from mfobs.checks.chlorination import _molecule
    for n, allowed in ((0, False), (1, True)):
        molecule = _molecule(0, n=n)
        dg = mf.DG(graphDatabase=[molecule])
        with dg.build() as builder:
            matches = builder.apply([molecule], rule, mode=mf.ApplyMode.UseAny, verbosity=0)
        assert bool(matches) is allowed
    with pytest.raises(UnsupportedCondition, match="Nonlocal"):
        LocalConditionalReducer().reduce(rule)


@pytest.mark.parametrize("condition", [
    'not [ exists [ node [ id 2 label "C" ] edge [ source 0 target 2 label ":" ] ] ]',
    'or [ exists [ node [ id 2 label "C" ] edge [ source 0 target 2 label ":" ] ] '
    'not [ exists [ node [ id 2 label "Cl" ] edge [ source 0 target 2 label "-" ] ] ] ]',
])
def test_ring_typing_cannot_be_inferred_from_negation_or_one_or_branch(condition):
    import molfoundry as mf
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.models.chlorination.counting import LocalConditionalReducer
    rule = mf.Rule.fromGMLString('rule [ context [ node [ id 0 label "C" ] '
        'node [ id 1 label "H" ] edge [ source 0 target 1 label "-" ] ] '
        'condition [ ' + condition + ' ] ]', add=False)
    with pytest.raises(UnsupportedCondition, match="Nonlocal"):
        LocalConditionalReducer().reduce(rule)


def test_growth_certificate_rejects_a_ring_changing_rule():
    from mfobs.core.conditions import UnsupportedCondition
    from mfobs.models.chlorination.events import _certify_growth
    backend = Chlorination()
    with pytest.raises(UnsupportedCondition, match="Growth may change"):
        _certify_growth(backend.rule, backend.patterns)
