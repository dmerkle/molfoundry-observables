"""Actual synthesis and a deliberately honest positive-encoding baseline."""
from mfobs.core.graph import edge_rows, graph_key, make_graph, node_labels, node_rows
from mfobs.core.rules import make_rule, rule_graph
from mfobs.core.gml import bond_order_sums
from mfobs.core.conditions import engine_conditions, engine_no_edges
from pathlib import Path

import pytest

from mfobs.core.conditions import matches, transitions

from mfobs.core.worklist import run
from mfobs.inputs import load
from mfobs.models import tethered_chemistry as chemistry

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('variant,dimension', [('guarded', 2), ('ablated', 1)])
def test_real_tethered_compositions_close_from_one_seed(variant, dimension):
    result = run(load(ROOT/f'examples/tethered-nac/model-{variant}.json'), max_expansions=4)
    assert result['closed'] and not result['unknown'] and not result['pending']
    assert result['initial_observables'] == ['B']
    assert len(result['basis']) == dimension
    assert result['preservation_certificate']['status'] == 'VERIFIED_ALL_MULTISETS'
    assert result['ordered_composition_calls'] == 2*dimension
    if variant == 'guarded':
        assert result['basis_counting_functions'] == {'B': {'A': '4', 'B': '4'}, 'O1': {'A': '2'}}
        assert result['matrix'] == [['0', '-2*k'], ['0', '-k']]
        assert result['basis_origins']['O1']['kind'] == 'computed-channel-residual'
    else:
        assert result['matrix'] == [['-k']]


def test_positive_hydrogen_context_replaces_guard_on_declared_family():
    # This comparison does NOT modify the actual conditional synthesis input.
    left, right = chemistry.reaction_sides()
    def extend(graph):
        return make_graph({**node_labels(graph), 6: 'H'}, [*edge_rows(graph), (0, 6, '-')])
    positive = make_rule('positive-comparison-only', extend(left), extend(right))
    guarded = chemistry.guarded_rule()
    assert set(engine_no_edges(positive)) == set(chemistry.PRODUCT_NO_EDGES)
    assert engine_conditions(positive) == ()
    assert node_rows(rule_graph(positive, "left")) == node_rows(rule_graph(positive, "right"))
    assert bond_order_sums(rule_graph(positive, "left")) == bond_order_sums(rule_graph(positive, "right"))
    populations = [dict(A=a, B=b, P=1, D=1) for a, b in [(0, 1), (1, 0), (1, 2), (2, 2)]]
    for population in populations:
        host = chemistry.host(population)
        guarded_maps = matches(guarded, host)
        positive_maps = matches(rule_graph(positive, "left"), host)
        # Restriction to the original six atoms is a bijection of matches.
        assert {tuple(sorted(mt.items())) for mt in guarded_maps} == {
            tuple(sorted((u, v) for u, v in mt.items() if u != 6)) for mt in positive_maps}
        assert len(guarded_maps) == len(positive_maps) == 2*population['A']
        assert {graph_key(g): n for g, n in transitions(guarded, host).items()} == {
            graph_key(g): n for g, n in transitions(positive, host).items()}
