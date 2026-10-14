"""Exact all-inventory linear counting and automatic connected-jump discovery."""
from mfobs.core.graph import node_rows, make_graph, edge_rows
from mfobs.core.rules import make_rule, rule_graph
import json

import pytest

from mfobs.core.conditions import Adjacency, engine_no_edges
from mfobs.core.gml import rule_gml

from mfobs.core.worklist import run
from mfobs.inputs import InputError, load
from mfobs.models.finite_components import Counter, ONE, identity


def graph_data(graph):
    return {'nodes': list(node_rows(graph)), 'edges': list(edge_rows(graph))}


def fixture_graphs():
    nodes = {0: 'A', 1: 'B', 2: 'C'}
    active = make_graph(nodes, [(0, 1, '='), (1, 2, '-')])
    blocked = make_graph(nodes, [*edge_rows(active), (0, 2, '-')])
    product = make_graph(nodes, [(0, 1, '-'), (1, 2, '-')])
    return {'active': active, 'blocked': blocked, 'product': product}


def language(graphs=None):
    return {'kind': 'finite-component-multiset', 'components': {
        name: graph_data(graph) for name, graph in (graphs or fixture_graphs()).items()}}


def reaction():
    graphs = fixture_graphs()
    return make_rule('convert', graphs['active'], graphs['product'], ((0, 2),))


def test_connected_count_is_exact_linear_form_with_context_guard():
    counter = Counter(language())
    seed = make_graph({0: 'A', 1: 'B'}, [(0, 1, '=')])
    assert counter.count(seed)['polynomial'] == {'active': 1, 'blocked': 1}
    rule = reaction()
    assert counter.count(rule_graph(rule), engine_no_edges(rule))['polynomial'] == {'active': 1}
    degree = Adjacency(0, ('B', 'C'), ('-', '='), '=', 1)
    assert counter.count(seed, adjacency=(degree,))['polynomial'] == {'active': 1}
    assert counter.count(make_graph({}))['polynomial'] == {ONE: 1}


def test_disconnected_nonzero_pattern_fails_closed():
    counter = Counter(language())
    disconnected = make_graph({0: 'A', 1: 'B'})
    assert counter.count(disconnected)['status'] == 'unknown'
    contradiction = make_graph({0: 'A', 1: 'B', 2: 'C'}, [(0, 1, '=')])
    assert counter.count(contradiction, ((0, 1),))['status'] == 'zero'


@pytest.mark.parametrize('invalid', ['empty', 'disconnected', 'duplicate'])
def test_component_library_rejects_ambiguous_or_noncomponent_hosts(invalid):
    graphs = fixture_graphs()
    if invalid == 'empty':
        graphs['bad'] = make_graph({})
    elif invalid == 'disconnected':
        graphs['bad'] = make_graph({0: 'A', 1: 'B'})
    else:
        graphs['bad'] = make_graph({u+10: label for u, label in node_rows(graphs['active'])},
                                   [(u+10, v+10, b) for u, v, b in edge_rows(graphs['active'])])
    with pytest.raises(InputError):
        Counter(language(graphs))


def test_preservation_checks_every_species_not_a_population_cutoff():
    certificate = Counter(language()).preservation({'convert': reaction()})
    assert certificate['status'] == 'VERIFIED_ALL_MULTISETS'
    assert {row['component']: row['admissible_applications']
            for row in certificate['checks']} == {'active': 1, 'blocked': 0, 'product': 0}
    incomplete = fixture_graphs(); del incomplete['product']
    with pytest.raises(InputError, match='not preserved'):
        Counter(language(incomplete)).preservation({'convert': reaction()})


def test_disconnected_rule_and_vertex_creation_are_rejected():
    counter = Counter(language())
    left = make_graph({0: 'A', 1: 'B'})
    with pytest.raises(InputError, match='nonempty and connected'):
        counter.preservation({'bad': identity(left)})
    left = make_graph({0: 'A'})
    right = make_graph({0: 'A', 1: 'B'}, [(0, 1, '-')])
    with pytest.raises(InputError, match='preserve every vertex'):
        counter.preservation({'bad': make_rule('bad', left, right)})


def test_composition_discovers_two_functions_without_target_basis(tmp_path):
    seed = make_graph({0: 'A', 1: 'B'}, [(0, 1, '=')])
    for name, value in {'rule.gml': rule_gml(reaction()),
                        'seed.gml': identity(seed).getGMLString()}.items():
        (tmp_path/name).write_text(value, encoding='utf-8')
    (tmp_path/'language.json').write_text(json.dumps(language()), encoding='utf-8')
    (tmp_path/'reduction.json').write_text(json.dumps({
        'method': 'exact-connected-components-and-dynamic-span'}), encoding='utf-8')
    manifest = {'schema_version': 1, 'id': 'test', 'title': 'Connected baseline',
        'adapter': 'mfobs.models.finite_components:Adapter',
        'rules': [{'id': 'convert', 'file': 'rule.gml', 'rate': 'k', 'per_match': '1', 'kind': 'abstract'}],
        'observables': [{'id': 'B', 'file': 'seed.gml'}], 'initial_observables': ['B'],
        'composition': {'mode': 'public', 'connected': False, 'include_empty': True},
        'forbidden': None, 'state_language': 'language.json', 'reduction': 'reduction.json',
        'normalisation': {}, 'initial_state': {'populations': {'active': 3, 'blocked': 2, 'product': 0}},
        'parameters': None, 'reference': None}
    (tmp_path/'model.json').write_text(json.dumps(manifest), encoding='utf-8')
    result = run(load(tmp_path/'model.json'), max_expansions=5)
    assert result['closed'] and result['basis'] == ['B', 'O1']
    assert result['basis_counting_functions'] == {
        'B': {'active': '1', 'blocked': '1'}, 'O1': {'active': '1'}}
    assert result['matrix'] == [['0', '-k'], ['0', '-k']]
    assert result['initial_values'] == {'B': '5', 'O1': '3'}
    assert result['basis_origins']['O1']['kind'] == 'computed-channel-residual'
    assert result['preservation_certificate']['status'] == 'VERIFIED_ALL_MULTISETS'
