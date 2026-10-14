"""Independent chemistry and direct-application checks for the tethered case."""
from mfobs.core.molecules import formula as molecular_formula
from mfobs.core.gml import bond_order_sums
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import edge_rows, has_edge, node_labels, node_rows
from mfobs.core.rules import rule_graph
from collections import Counter

import networkx as nx
import pytest

from mfobs.core.conditions import matches, transitions
from mfobs.models import tethered_chemistry as chemistry


def isomorphic(a, b):
    def convert(graph):
        result = nx.Graph()
        result.add_nodes_from((u, {'label': label}) for u, label in node_rows(graph))
        result.add_edges_from((u, v, {'label': label}) for u, v, label in edge_rows(graph))
        return result
    return nx.is_isomorphic(convert(a), convert(b),
                            node_match=nx.algorithms.isomorphism.categorical_node_match('label', None),
                            edge_match=nx.algorithms.isomorphism.categorical_edge_match('label', None))


def test_all_four_hosts_are_ordinary_covalent_and_saturated():
    hosts = chemistry.language_hosts()
    assert set(hosts) == {'A', 'B', 'P', 'D'}
    for name, graph in hosts.items():
        hydrogen = 10 if name in ('A', 'P') else 8
        assert molecular_formula(graph) == {'C': 6, 'H': hydrogen}
        assert graph.numVertices == 6 + hydrogen
        assert set(node_labels(graph).values()) == {'C', 'H'}
        assert {label for _, _, label in edge_rows(graph)} <= {'-', '='}
        sums = bond_order_sums(graph)
        assert all(sums[u] == (4 if label == 'C' else 1) for u, label in node_rows(graph))
        mol = to_rdkit(graph)
        assert mol.GetNumAtoms() == graph.numVertices
        assert all(atom.GetFormalCharge() == 0 and atom.GetNumRadicalElectrons() == 0
                   for atom in mol.GetAtoms())


def test_rule_atom_identity_valence_and_active_context_guard():
    rule, ablated = chemistry.guarded_rule(), chemistry.ablated_rule()
    assert node_rows(rule_graph(rule, "left")) == node_rows(rule_graph(rule, "right"))
    assert bond_order_sums(rule_graph(rule, "left")) == bond_order_sums(rule_graph(rule, "right"))
    for name, graph in chemistry.language_hosts().items():
        raw = len(matches(rule_graph(rule, "left"), graph))
        assert raw == (2 if name in ('A', 'B') else 0)
        guarded_targets = transitions(rule, graph)
        ablated_targets = transitions(ablated, graph)
        assert sum(guarded_targets.values()) == (2 if name == 'A' else 0)
        assert sum(ablated_targets.values()) == (2 if name in ('A', 'B') else 0)
        for target in (*guarded_targets, *ablated_targets):
            assert node_rows(target) == node_rows(graph)
            assert bond_order_sums(target) == bond_order_sums(graph)
            expected = chemistry.product() if name == 'A' else chemistry.blocked_product()
            assert isomorphic(target, expected)
    assert has_edge(chemistry.blocked(), 0, 2, '-')
    assert has_edge(chemistry.blocked_product(), 0, 2, '-')


@pytest.mark.parametrize('population', [
    {}, {'A': 1}, {'B': 1}, {'A': 1, 'B': 1},
    {'A': 2, 'B': 3, 'P': 2, 'D': 1}, {'B': 2, 'P': 1},
])
def test_arbitrary_mixture_counts_and_seed(population):
    graph = chemistry.host(population)
    a, b = population.get('A', 0), population.get('B', 0)
    assert len(matches(chemistry.seed_graph(), graph)) == 4*(a+b)
    assert sum(transitions(chemistry.guarded_rule(), graph).values()) == 2*a
    assert sum(transitions(chemistry.ablated_rule(), graph).values()) == 2*(a+b)
    assert chemistry.PER_MATCH * sum(transitions(chemistry.guarded_rule(), graph).values()) == a


def test_no_cross_molecule_reaction_and_closed_products():
    # The connected positive tether cannot be assembled across separate hosts.
    for graph in (chemistry.product(), chemistry.blocked_product(),
                  chemistry.host({'P': 3, 'D': 2})):
        assert not transitions(chemistry.guarded_rule(), graph)
        assert not transitions(chemistry.ablated_rule(), graph)
    with pytest.raises(ValueError):
        chemistry.host({'unknown': 1})
    with pytest.raises(ValueError):
        chemistry.host({'A': -1})


def test_public_molfoundry_targets_with_explicit_finite_inventory():
    from mfobs.core.gml import component_graphs, graph_gml, rule_gml
    import molfoundry as mf
    def keys(graph):
        return tuple(sorted(mf.Graph.fromGMLString(graph_gml(part), add=False).canonKey
                            for part in component_graphs(graph)))
    checked = 0
    observed = []
    for population in ({'A': 1}, {'B': 1}, {'P': 1}, {'D': 1},
                       {'A': 2, 'B': 1, 'P': 1}):
        graph = chemistry.host(population)
        parts = [mf.Graph.fromGMLString(graph_gml(part), add=False)
                 for part in component_graphs(graph)]
        inventory = Counter(part.canonKey for part in parts)
        species = list({part.canonKey: part for part in parts}.values())
        for rule in (chemistry.guarded_rule(), chemistry.ablated_rule()):
            engine_rule = mf.Rule.fromGMLString(rule_gml(rule), add=False)
            dg = mf.DG(graphDatabase=species)
            actual = set()
            with dg.build() as builder:
                for edge in builder.apply(species, engine_rule, mode=mf.ApplyMode.UseAnyMultiplicity):
                    consumed = Counter(v.graph.canonKey for v in edge.sources)
                    if any(count > inventory[key] for key, count in consumed.items()):
                        continue
                    remaining = inventory.copy()
                    remaining.subtract(consumed)
                    remaining.update(v.graph.canonKey for v in edge.targets)
                    actual.add(tuple(sorted(remaining.elements())))
            observed.append((population, actual, rule, graph))
            checked += 1
    assert checked == 10
    for label, actual, rule, graph in observed:
        assert actual == {keys(target) for target in transitions(rule, graph)}, label
