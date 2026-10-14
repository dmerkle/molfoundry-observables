"""Applications at a given match through MolFoundry's public applyMatch."""
from mfobs.core.graph import same_graph, graph_key, make_graph, edge_rows, node_labels, apply_match, matches, transitions
from mfobs.core.rules import make_rule, rule_graph
from collections import Counter

import molfoundry as mf
import pytest


# A methyl group and a Y on an X.
HOST = make_graph({10: 'X', 11: 'C', 12: 'H', 13: 'H', 14: 'H', 20: 'Y'},
                  [(10, 11, '-'), (11, 12, '-'), (11, 13, '-'), (11, 14, '-'), (10, 20, '-')])
X = make_graph({0: 'X'})
METHYL_LEAVES = make_rule('methyl-leaves',
                     make_graph({0: 'X', 1: 'C', 2: 'H', 3: 'H', 4: 'H'},
                                [(0, 1, '-'), (1, 2, '-'), (1, 3, '-'), (1, 4, '-')]),
                     X)
ADD_Z = make_rule('add-z', X, make_graph({0: 'X', 5: 'Z'}, [(0, 5, '-')]))


def test_product_keeps_host_ids_and_numbers_created_vertices():
    product = apply_match(ADD_Z, HOST, {0: 10})
    # Z takes the smallest ID the host does not use.
    assert same_graph(product, make_graph({**node_labels(HOST), 0: 'Z'}, [*edge_rows(HOST), (0, 10, '-')]))


def test_inadmissible_matches_have_no_product():
    # Deleting X alone would leave its edges dangling.
    assert apply_match(make_rule('delete-x', X, make_graph({})), HOST, {0: 10}) is None
    assert apply_match(ADD_Z, HOST, {0: 99}) is None
    assert apply_match(ADD_Z, HOST, {0: 11}) is None  # a C is not an X


def test_editing_an_enumerated_match_requires_native_validation():
    host = make_graph({0: 'C', 1: 'O'}, [(0, 1, '-')])
    rule = make_rule('delete-pair', host, make_graph({}))
    match = matches(rule, host)[0]
    assert same_graph(apply_match(rule, host, match), make_graph({}))
    match[0], match[1] = match[1], match[0]
    assert apply_match(rule, host, match) is None


def test_engine_failure_on_an_enumerated_rule_match_propagates(monkeypatch):
    match = matches(ADD_Z, HOST)[0]

    def fail(*args):
        raise mf.LogicError('unexpected engine failure')

    monkeypatch.setattr(mf.DG.Builder, 'applyMatch', fail)
    with pytest.raises(mf.LogicError, match='unexpected engine failure'):
        apply_match(ADD_Z, HOST, match)


def test_left_pattern_match_is_not_a_validated_rule_application():
    delete = make_rule('delete-x', X, make_graph({}))
    assert apply_match(delete, HOST, matches(X, HOST)[0]) is None


def test_transitions_count_every_match_of_a_shared_application():
    assert len(matches(rule_graph(METHYL_LEAVES), HOST)) == 6  # the hydrogens permute
    leaves = make_graph({10: 'X', 20: 'Y'}, [(10, 20, '-')])
    assert Counter({graph_key(g): n for g, n in transitions(METHYL_LEAVES, HOST).items()}) == Counter({graph_key(leaves): 6})


@pytest.mark.parametrize('host', [
    make_graph({}),
    make_graph({0: 'A', 2: 'B'}, [(0, 2, '-')]),
])
def test_birth_uses_public_application_and_preserves_host(host, monkeypatch):
    # Two isomorphic products must stay distinct, even from an identical host.
    right = make_graph({5: 'A', 6: 'B', 7: 'A', 8: 'B'},
                       [(5, 6, '-'), (7, 8, '-')])
    birth = make_rule('birth', make_graph({}), right)
    applications = []
    original = mf.DG.Builder.applyMatch

    def record_application(builder, graphs, rule, match):
        applications.append(dict(match))
        return original(builder, graphs, rule, match)

    monkeypatch.setattr(mf.DG.Builder, 'applyMatch', record_application)
    product = apply_match(birth, host, {})
    fresh = [node for node in range(6) if node not in node_labels(host)][:4]
    assert same_graph(product, make_graph(
        {**node_labels(host), **dict(zip(fresh, ['A', 'B', 'A', 'B']))},
        [*edge_rows(host), (fresh[0], fresh[1], '-'), (fresh[2], fresh[3], '-')],
    ))
    assert applications == [{}]


@pytest.mark.parametrize('host', [make_graph({}), HOST])
def test_empty_rule_preserves_host(host):
    empty = make_graph({})
    assert same_graph(apply_match(make_rule('empty', empty, empty), host, {}), host)
