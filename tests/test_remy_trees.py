
from mfobs.core.graph import edge_rows, node_labels, make_graph, matches
import mfobs.models.remy.analysis as remy_analysis
import mfobs.models.remy.counting as remy_counting
import mfobs.models.remy.events as remy_events
import mfobs.checks.remy_trees as remy_remy_trees
from mfobs.checks import remy_trees as remy_check
from mfobs.numerics import remy_trees as remy_numerics
"""Independent and regression checks for the Rémy-tree artifact."""
from collections import Counter
from fractions import Fraction
from math import factorial
from pathlib import Path
import json

import sympy as sp

from mfobs.core.engine import mf, sqpo_jump_closure
from mfobs.core.gml import rule_side
from mfobs.core.rules import rule_graph, rule_vertex_ids, removed_edges, is_identity
from mfobs.models.remy import inputs as remy

ROOT = Path(__file__).resolve().parents[1]


def _plain_tree(tree):
    nodes, typed_edges = remy._typed_tree_graph(tree)
    edges = []
    edge_id = len(nodes)
    for source, target, edge_type in typed_edges:
        nodes.append((edge_id, edge_type))
        edges.extend([(source, edge_id, "src"), (edge_id, target, "tgt")])
        edge_id += 1
    return make_graph(dict(nodes), edges)


def test_six_typed_rules_preserve_direction_order_and_endpoints():
    rules = remy.insertion_rules()
    assert list(rules) == [
        "GL_I", "GL_L", "GL_R", "GR_I", "GR_L", "GR_R"
    ]
    assert len({rule.canonKey for rule in rules.values()}) == 6
    for name, rule in rules.items():
        assert set(node_labels(rule_graph(rule, "context"))) == {0, 1}
        assert len(rule_vertex_ids(rule, "deleted")) == 1
        assert len(rule_vertex_ids(rule, "created")) == 5
        assert {label for _, _, label in removed_edges(rule)} == {"src", "tgt"}
        right = node_labels(rule_graph(rule, "right"))
        right_labels = Counter(right[node] for node in rule_vertex_ids(rule, "created"))
        assert right_labels[name[-1]] >= 1
        assert right_labels["L"] >= 1 and right_labels["R"] >= 1


def test_molfoundry_gml_edges_are_natively_undirected():
    forward = '''graph [
      node [ id 0 label "A" ] node [ id 1 label "B" ]
      edge [ source 0 target 1 label "x" ]
    ]'''
    reverse = '''graph [
      node [ id 0 label "A" ] node [ id 1 label "B" ]
      edge [ source 1 target 0 label "x" ]
    ]'''
    assert (
        mf.graphGMLString(forward, add=False).canonKey
        == mf.graphGMLString(reverse, add=False).canonKey
    )


def test_observable_gml_counts_recursive_patterns_on_actual_incidence_graphs():
    rules = remy.observable_rules()
    for size in range(5):
        for tree in remy_check.trees(size):
            host = _plain_tree(tree)
            actual = []
            for name in remy.OBSERVABLES:
                count = 0
                for component in rules[name]:
                    pattern = rule_side(component, "left")
                    count += len(matches(pattern, host))
                actual.append(count)
            assert tuple(actual) == remy_check.observable_values(tree)


def test_symbolic_all_tree_derivation_reproduces_eq57g_without_enumeration():
    result = remy_events.derive_jump_classes()
    assert result["status"] == "VERIFIED_ALL_VALID_TREES"
    assert not result["enumeration_used"]
    actual = {
        tuple(event["delta"]): tuple(event["propensity"])
        for event in result["classes"]
    }
    # Published Eq. 57g is a validation oracle, never an analyzer input.
    assert actual == {
        (2, 1, 0, 0): (2, -3, 2, -1),
        (2, 1, 1, 0): (0, 3, -6, 3),
        (2, 1, 1, 1): (0, 0, 4, -3),
        (2, 1, 0, -1): (0, 0, 0, 1),
    }
    assert all(
        len(event["propensity"]) == 4 for event in result["classes"]
    )


def test_generic_worklist_discovers_four_dimensional_span_and_auxiliary_identity():
    analysis = remy_analysis.operator_analysis()
    span = analysis["span"]
    assert span["status"] == "CLOSED" and span["dimension"] == 4
    vectors = {
        row["expression"]: tuple(map(int, row["event_coefficients"]))
        for row in span["basis"]
    }
    generator = vectors["G"]
    p2g = vectors["ad_P2(G)"]
    p3g = vectors["ad_P3(G)"]
    p3p2g = vectors["ad_P3(ad_P2(G))"]
    auxiliary = tuple(a - b for a, b in zip(p3p2g, p3g))
    jumps = {
        name: tuple(event["delta"][index] for event in remy_events.derive_jump_classes()["classes"])
        for index, name in enumerate(remy.OBSERVABLES)
    }
    action = lambda observable, value: tuple(
        jump * coefficient for jump, coefficient in zip(jumps[observable], value)
    )
    assert action("E", generator) == tuple(2 * value for value in generator)
    assert action("P1", generator) == generator
    assert action("P2", p2g) == p2g
    assert action("P2", auxiliary) == (0, 0, 0, 0)
    assert action("P3", auxiliary) == tuple(-value for value in auxiliary)


def test_generated_emgf_is_eq57g_and_eq58_monomials_follow_jumps():
    classes = remy_events.derive_jump_classes()["classes"]
    epsilon, gamma, mu, nu = sp.symbols("epsilon gamma mu nu")
    variables = (epsilon, gamma, mu, nu)
    common = sp.exp(2 * epsilon + gamma)
    coefficients = []
    for observable in range(4):
        value = sum(
            event["propensity"][observable]
            * sp.exp(sum(d * variable for d, variable in zip(event["delta"], variables)))
            for event in classes
        )
        coefficients.append(sp.simplify(value / common))
    expected = (
        2,
        3 * (sp.exp(mu) - 1),
        4 * sp.exp(mu + nu) - 6 * sp.exp(mu) + 2,
        3 * sp.exp(mu) + sp.exp(-nu) - 3 * sp.exp(mu + nu) - 1,
    )
    assert all(sp.simplify(a - b) == 0 for a, b in zip(coefficients, expected))
    monomials = [
        value["factor"]
        for value in remy_analysis.operator_analysis()["variable_substitution_monomials"]
    ]
    assert monomials == [
        "x_E^2*x_P1",
        "x_E^2*x_P1*x_P2",
        "x_E^2*x_P1*x_P2*x_P3",
        "x_E^2*x_P1/x_P3",
    ]


def test_finite_oracle_and_exact_joint_distributions_are_validation_only():
    validation = remy_check.validate_finite()
    assert validation["source_trees"] == 2056
    assert validation["rewrite_events_counting_multiplicity"] == 66196
    distributions = remy_numerics.joint_distributions((3, 4, 5))
    for generation in (3, 4, 5):
        expected = Counter(
            remy_check.observable_values(tree)[2:]
            for tree in remy_check.trees(generation)
        )
        actual = {
            (row["P2"], row["P3"]): int(row["weight"])
            for row in distributions["generations"][str(generation)]["states"]
        }
        per_tree_weight = factorial(generation + 1)
        assert actual == {
            key: count * per_tree_weight for key, count in expected.items()
        }


def test_existing_assessment_oracle_agrees_but_is_not_a_synthesis_input():
    oracle = json.loads(
        (ROOT / "examples/remy-trees/finite-oracle.json").read_text(encoding="utf-8")
    )
    actual = {
        tuple(event["delta"]): tuple(map(str, event["propensity"]))
        for event in remy_events.derive_jump_classes()["classes"]
    }
    expected = {
        tuple(event["delta"]): tuple(
            event["propensity_coefficients_E_P1_P2_P3"]
        )
        for event in oracle["jump_classes"]
    }
    assert actual == expected
    assert "validation only" in oracle["provenance"]


def test_actual_molfoundry_rules_have_the_same_direct_targets_and_match_counts():
    result = remy_check.validate_molfoundry_direct(3)
    assert result["status"] == "PASS"
    assert result["injective_LHS_matches"] == 98
    assert result["match_multiplicities_preserved"]


def test_sqpo_jump_closure_is_identity_on_full_input():
    rule = remy.insertion_rules()["GL_I"]
    jump = sqpo_jump_closure(rule)
    assert is_identity(jump)
    kept = rule_graph(jump, "context")
    assert len(node_labels(kept)) == 3
    assert len(edge_rows(kept)) == 2


def test_shared_composition_compiler_derives_the_constant_remy_prefix():
    result = remy_events.composition_prefix_compilation()
    assert result.evidence["jump_projection"] == remy_counting.REMY_JUMP_PROJECTION
    assert result.model.provenance.all_jumps_composition_derived
    assert [(event.identifier, event.delta) for event in result.model.events] == [
        ("GL_I", (2, 1)),
        ("GL_L", (2, 1)),
        ("GL_R", (2, 1)),
        ("GR_I", (2, 1)),
        ("GR_L", (2, 1)),
        ("GR_R", (2, 1)),
    ]
    assert [result.model.ring.text(event.propensity) for event in result.model.events] == [
        "1", "P1", "P1", "1", "P1", "P1",
    ]


def test_observables_are_synthesized_from_p3_alone(monkeypatch):
    original = remy.observable_gml
    requested = []

    def recording_gml(name, stem):
        requested.append((name, stem))
        return original(name, stem)

    monkeypatch.setattr(remy_counting, "observable_gml", recording_gml)
    result = remy_counting.observable_synthesis()
    assert result.status == "CLOSED"
    assert list(result.basis) == ["P3", "P2", "P1", "E"]
    assert requested == [("P3", stem) for stem in remy.EDGE_TYPES]
    assert [
        (step["source"], step["coordinates"], step["added"])
        for step in result.steps
    ] == [
        ("P3", {"P3": "-4", "P2": "4"}, "P2"),
        ("P2", {"P2": "-2", "P1": "3"}, "P1"),
        ("P1", {"E": "2"}, "E"),
        ("E", {"E": "4"}, None),
    ]
    assert all(
        origin["origin"] == "computed-channel-residual"
        for name, origin in result.origins.items()
        if name != "P3"
    )


def test_raw_molfoundry_commutators_retain_terms_and_overlap_rejections():
    raw = remy_remy_trees.raw_first_commutators()
    assert {name: len(value["coefficients"]) for name, value in raw.items()} == {
        "E": 6, "P1": 6, "P2": 15, "P3": 24,
    }
    assert set(raw["E"]["coefficients"].values()) == {2}
    assert set(raw["P1"]["coefficients"].values()) == {1}
    assert sum(
        row[side]["forbidden_overlap_rejected"]
        for value in raw.values()
        for row in value["ordered_composition_rows"]
        for side in ("left", "right")
    ) > 0
    reasons = Counter(
        reason
        for value in raw.values()
        for row in value["ordered_composition_rows"]
        for side in ("left", "right")
        for reason, count in row[side]["forbidden_reasons"].items()
        for _ in range(count)
    )
    assert sum(
        count for reason, count in reasons.items()
        if reason.endswith(":directed-cycle")
    ) > 0


def test_guarded_dpo_reduction_covers_discovery_and_cancelled_composites():
    """Audit the rewrite terms, not just the six input rules or their totals."""
    from mfobs.core.engine import dpo_jump_closure
    from mfobs.core.identity import jump_input
    from mfobs.core.typed_incidence import typed_incidence_commutator

    representatives = dict((v.canonKey,v) for v in remy.insertion_rules().values())
    synthesis = remy_counting.observable_synthesis()
    assert synthesis.jump_projection == remy_counting.REMY_JUMP_PROJECTION
    assert list(synthesis.basis) == ['P3','P2','P1','E']
    for evidence in synthesis.evidence:
        representatives.update(evidence['raw'].representatives)
    for components in remy.observable_rules().values():
        for component in components:
            for reaction in remy.insertion_rules().values():
                pair = typed_incidence_commutator(component,reaction,
                    specification=remy.TYPED_INCIDENCE_SPECIFICATION,
                    connected=False,include_empty=True)
                # Includes representatives whose coefficients cancel to zero.
                representatives.update(pair.representatives)
    assert len(representatives) == 384
    for rule in representatives.values():
        projected = dpo_jump_closure(rule)
        deleted = rule_vertex_ids(projected, "deleted")
        # I \ K is deleted and recreated under a fresh id (Thm. 6(2)).
        assert len(deleted) == 1 and len(rule_vertex_ids(projected, "created")) == 1
        assert jump_input(projected)[1] == deleted
        reduced = remy_counting.valid_tree_dpo_jump_closure(rule)
        assert reduced.canonKey == sqpo_jump_closure(rule).canonKey
        assert is_identity(reduced)
    assert 'DPO' in remy_events.polynomial_jump_model().provenance.jump_projection
    assert 'DPO' in remy_analysis.operator_analysis()['emgf_operator']['construction']


def test_valid_tree_reduction_rejects_unsupported_deletions_and_conditions():
    import pytest
    bad = [
        ('rule [ left [ node [ id 0 label "V" ] ] ]','complete src/tgt'),
        ('rule [ context [ node [ id 0 label "V" ] ] left [ node [ id 1 label "I" ] edge [ source 0 target 1 label "src" ] ] ]','complete src/tgt'),
        ('rule [ context [ node [ id 0 label "V" ] node [ id 1 label "I" ] ] left [ edge [ source 0 target 1 label "src" ] ] ]','complete edge-nodes'),
        ('rule [ context [ node [ id 0 label "V" ] node [ id 1 label "V" ] ] constrainNoEdge [ source 0 target 1 ] ]','additional conditions'),
    ]
    for gml,message in bad:
        with pytest.raises(ValueError,match=message):
            remy_counting.valid_tree_dpo_jump_closure(mf.Rule.fromGMLString(gml,add=False))


def test_linear_discovery_does_not_need_the_pk_recognizer(monkeypatch):
    """Disable structural recognition, rather than just relabeling its output."""
    from mfobs.core.jump_closure_worklist import ObservableDesignation

    named = remy_counting.observable_synthesis(order=3)
    monkeypatch.setattr(
        remy_counting, 'classify_observable',
        lambda expression, index: ObservableDesignation(f'O{index}'),
    )
    neutral = remy_counting.observable_synthesis(order=3)
    assert neutral.status == named.status == 'CLOSED'
    assert list(neutral.basis) == ['P3', 'O1', 'O2', 'O3']
    assert [e.coefficients for e in neutral.basis.values()] == [
        e.coefficients for e in named.basis.values()
    ]
    names = dict(zip(neutral.basis, named.basis))
    for before, after in zip(named.steps, neutral.steps, strict=True):
        assert names[after['source']] == before['source']
        assert {names[k]: v for k, v in after['coordinates'].items()} == before['coordinates']
        assert names.get(after['added']) == before['added']


def test_p4_linear_closure_survives_incomplete_event_reduction():
    from mfobs import RemyTrees, analyze
    backend = RemyTrees(order=4)
    result = analyze(backend.specification(), backend=backend, linear_limit=8).record()
    assert result['linear']['status'] == 'CLOSED'
    assert result['linear']['basis'] == ['P4', 'P3', 'P2', 'P1', 'E']
    assert result['pjc']['status'] == 'UNKNOWN_UNSUPPORTED_EVENT_REDUCTION'
    assert any('P5' in row['required_observables'] for row in result['pjc']['unresolved'])
    assert 'operator' not in result['pjc']
    assert 'classes' not in remy_events.derive_jump_classes(4)

    # Count the query definition on concrete trees, then apply every insertion.
    # This does not use the context partition or its local jump formula.
    def values(tree):
        if not tree:
            return {'E': 1, **{f'P{k}': 0 for k in range(1, 5)}}
        left, right = map(values, tree)
        depth, branch = 0, tree
        while branch:
            depth += 1
            branch = branch[0]
        counts = {f'P{k}': left[f'P{k}'] + right[f'P{k}'] + (depth >= k)
                  for k in range(1, 5)}
        return {'E': 2 * counts['P1'] + 1, **counts}

    for size in range(7):
        for tree in remy_check.trees(size):
            before = values(tree)
            after = [values(event.target) for event in remy_check.insertion_events(tree)]
            for step in result['linear']['steps']:
                name = step['source']
                actual = sum(row[name] - before[name] for row in after)
                predicted = sum(Fraction(c) * before[n] for n, c in step['coordinates'].items())
                assert actual == predicted, (size, name, tree)
