"""Stable state IDs and native conditions across the public execution boundary."""
import molfoundry as mf

from mfobs.core.gml import graph_gml
from mfobs.core.graph import edge_rows, graph_key, make_graph, node_labels, same_graph
from mfobs.core.molfoundry_execution import apply_match, count_matches, matches, transitions
from mfobs.core.rules import make_rule, rule_graph


def test_graph_identity_is_distinct_from_exact_state_and_isomorphism_keys():
    first = make_graph({7: 'C', 20: 'O'}, [(7, 20, '=')])
    second = mf.Graph.fromGMLString(graph_gml(first), add=False)
    renamed = make_graph({0: 'C', 1: 'O'}, [(0, 1, '=')])
    assert isinstance(first, mf.Graph)
    assert first != second  # MolFoundry objects retain their own identity.
    assert same_graph(first, second)
    assert graph_key(first) != graph_key(renamed)
    assert first.canonKey == renamed.canonKey
    assert node_labels(second) == {7: 'C', 20: 'O'}
    assert edge_rows(second) == ((7, 20, '='),)


def test_native_nac_is_shared_by_count_enumeration_and_application():
    rule = mf.Rule.fromGMLString('''rule [
      context [ node [ id 7 label "C" ] ]
      right [ node [ id 8 label "H" ] edge [ source 7 target 8 label "-" ] ]
      condition [ not [ exists [ node [ id 9 label "O" ]
        edge [ source 7 target 9 label "-" ] ] ] ]
    ]''', add=False)
    allowed = make_graph({33: 'C'})
    forbidden = make_graph({33: 'C', 70: 'O'}, [(33, 70, '-')])
    assert count_matches(rule, allowed) == len(matches(rule, allowed)) == 1
    assert count_matches(rule, forbidden) == len(matches(rule, forbidden)) == 0
    assert apply_match(rule, forbidden, {7: 33}) is None
    product = apply_match(rule, allowed, {7: 33})
    assert isinstance(product, mf.Graph)
    assert node_labels(product) == {33: 'C', 0: 'H'}
    assert edge_rows(product) == ((0, 33, '-'),)


def test_native_rule_side_keeps_ids_and_symmetric_products_keep_multiplicity():
    rule = make_rule('bind', make_graph({7: 'A', 20: 'A'}),
                     make_graph({7: 'A', 20: 'A'}, [(7, 20, '-')]))
    assert isinstance(rule, mf.Rule)
    assert node_labels(rule_graph(rule)) == {7: 'A', 20: 'A'}
    host = make_graph({33: 'A', 70: 'A'})
    assert len(matches(rule, host)) == count_matches(rule, host) == 2
    products = transitions(rule, host)
    assert len(products) == 1
    product, multiplicity = next(iter(products.items()))
    assert multiplicity == 2
    assert graph_key(product) == graph_key(make_graph({33: 'A', 70: 'A'}, [(33, 70, '-')]))
