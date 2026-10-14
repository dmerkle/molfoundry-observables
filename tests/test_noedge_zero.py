"""Sound zero elimination without assumptions about the molecular language."""
from mfobs.core.graph import make_graph, node_labels
from mfobs.core.rules import make_rule
from pathlib import Path

from mfobs.core.conditions import Adjacency, matches, engine_no_edges
from mfobs.core.engine import mf
from mfobs.core.gml import rule_gml

from mfobs.core.identity import (
    identity_contradictions,
    noedge_contradictions,
    normalized_identity_rule,
)
from mfobs.core.worklist import run
from mfobs.inputs import load


def test_positive_edge_and_noedge_are_contradictory_for_any_bond_label():
    for bond in ['-', '=']:
        pattern = make_graph({0: 'C', 1: 'C'}, [(0, 1, bond)])
        rule = make_rule('contradiction', pattern, pattern, ((0, 1),))
        engine_rule = mf.Rule.fromGMLString(rule_gml(rule), add=False)
        assert noedge_contradictions(engine_rule) == ((0, 1),)
        for host in [pattern, make_graph(node_labels(pattern)),
                     make_graph({0: 'C', 1: 'C', 2: 'C'},
                                [(0, 1, bond), (1, 2, bond)])]:
            assert not matches(pattern, host, engine_no_edges(rule))


def test_distinct_absent_pair_is_not_discarded():
    pattern = make_graph({0: 'C', 1: 'C', 2: 'C'}, [(0, 1, '=')])
    rule = make_rule('valid-NAC', pattern, pattern, ((0, 2),))
    assert not noedge_contradictions(mf.Rule.fromGMLString(rule_gml(rule), add=False))
    assert matches(pattern, pattern, engine_no_edges(rule))


def test_product_only_edge_does_not_contradict_input_noedge():
    left = make_graph({0: 'C', 1: 'C'})
    right = make_graph(node_labels(left), [(0, 1, '-')])
    rule = make_rule('create-edge-with-safeguard', left, right, ((0, 1),))
    assert not noedge_contradictions(mf.Rule.fromGMLString(rule_gml(rule), add=False))


def test_adjacency_identity_solver_eliminates_only_an_empty_integer_interval():
    pattern = make_graph({0: "Position", 1: "Next"}, [(0, 1, "src")])
    positive = Adjacency(0, ("Next",), ("src",), ">", 0)
    zero = Adjacency(0, ("Next",), ("src",), "=", 0)
    contradictory = mf.Rule.fromGMLString(rule_gml(make_rule(
        "zero", pattern, pattern, adjacency=(positive, zero)
    )), add=False)
    reasons = identity_contradictions(contradictory)
    assert len(reasons) == 1
    assert reasons[0]["kind"] == "incompatible-adjacency-conditions"
    assert reasons[0]["lower_bound"] == 1
    assert reasons[0]["upper_bound"] == 0

    duplicate = mf.Rule.fromGMLString(rule_gml(make_rule(
        "duplicate", pattern, pattern, adjacency=(positive, positive)
    )), add=False)
    single = mf.Rule.fromGMLString(rule_gml(make_rule(
        "single", pattern, pattern, adjacency=(positive,)
    )), add=False)
    assert not identity_contradictions(duplicate)
    assert normalized_identity_rule(duplicate).canonKey == normalized_identity_rule(single).canonKey


def test_nonadjacent_worklist_records_zero_but_remains_unknown():
    root = Path(__file__).resolve().parents[1]
    result = run(load(root/'examples/nonadjacent-pairs/model.json'), max_expansions=1)
    assert result['status'] == 'UNKNOWN' and not result['closed']
    assert result['basis'] == ['B', 'O1']
    assert result['pending'] == ['O1']
    assert result['equations_by_reaction']['B']['cycloaddition'] == {'O1': '-1/2'}
    assert len(result['zero_jump_terms']) == 1
    assert next(iter(result['zero_jump_terms'].values()))['contradictory_pairs']
    assert result['ordered_composition_calls'] == 2
    assert result['raw_compositions'] == 26
