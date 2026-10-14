"""Molecular substrates and atom conservation in net ring opening."""
from mfobs.chemical_depiction import to_rdkit

import networkx as nx
import pytest

from mfobs.core.conditions import matches, transitions
from mfobs.core.graph import edge_rows, has_edge, make_graph, node_labels, node_rows
from mfobs.core.rules import make_rule
from mfobs.models.molecular_rop import inputs as chemistry
from mfobs.core.gml import bond_order_sums, disjoint_union
from mfobs.core.molecules import VALENCE, complete_with_h, formula


def isomorphic(left, right):
    def nx_graph(graph):
        result = nx.Graph()
        result.add_nodes_from((u, {'label': label}) for u, label in node_rows(graph))
        result.add_edges_from((u, v, {'label': label}) for u, v, label in edge_rows(graph))
        return result
    return nx.is_isomorphic(
        nx_graph(left), nx_graph(right),
        node_match=nx.algorithms.isomorphism.categorical_node_match('label', None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match('label', None))


def test_explicit_monomer_and_chain_formulas_and_valences():
    assert formula(chemistry.monomer("B")) == {'C': 6, 'H': 10, 'O': 2}
    assert formula(chemistry.chain(("B",) * 0)) == {'C': 7, 'H': 8, 'O': 1}
    for n in (0, 1, 2, 7, 30):
        graph = chemistry.chain(("B",) * n)
        assert formula(graph) == {'C': 7 + 6*n, 'H': 8 + 10*n, 'O': 1 + 2*n}
        assert set(label for _, _, label in edge_rows(graph)) == {'-', '='}
        sums = bond_order_sums(graph)
        assert all(sums[u] == VALENCE[label] for u, label in node_rows(graph))
        oh = make_graph({0: 'O', 1: 'H'}, [(0, 1, '-')])
        assert len(matches(oh, graph)) == 1
        # Monomer recognition contains the complete cyclic ester, not any ester.
        assert not matches(chemistry.lactone_core("B"), graph)
    with pytest.raises(ValueError, match='unknown repeat-unit symbol'):
        chemistry.chain_heavy(('C',))


def test_net_ring_opening_preserves_every_atom_and_its_bond_order_sum():
    left, right = chemistry.propagation_sides("B")
    assert node_rows(left) == node_rows(right)
    assert bond_order_sums(left) == bond_order_sums(right)
    assert has_edge(left, 8, 9, '-') and not has_edge(right, 8, 9)
    assert has_edge(right, 2, 9, '-') and not has_edge(left, 2, 9)
    assert has_edge(left, 0, 2, '-') and not has_edge(right, 0, 2)
    assert has_edge(right, 0, 8, '-') and not has_edge(left, 0, 8)


@pytest.mark.parametrize('n', [0, 1, 3])
@pytest.mark.parametrize('m', [0, 1, 2])
def test_independent_net_propagation_all_products_match_next_chain(n, m):
    rule = make_rule('test-net-ring-opening', *chemistry.propagation_sides("B"))
    host = disjoint_union([chemistry.chain(("B",) * n)] + [chemistry.monomer("B")]*m)
    targets = transitions(rule, host)
    assert sum(targets.values()) == m
    if m:
        expected = disjoint_union([chemistry.chain(("B",) * (n+1))] + [chemistry.monomer("B")]*(m-1))
        for target in targets:
            assert node_rows(target) == node_rows(host)
            assert bond_order_sums(target) == bond_order_sums(host)
            assert isomorphic(target, expected)


def test_cyclic_ester_pattern_rejects_acyclic_and_other_ring_size_esters():
    original = chemistry.monomer("B")
    # A bond-cut acyclic analogue contains an ester but not a full lactone ring.
    opened = make_graph(node_labels(original), [edge for edge in edge_rows(original) if edge != (0, 7, '-')])
    assert not matches(chemistry.lactone_core("B"), opened)
    # Delta-valerolactone is not epsilon-caprolactone; no wildcard polymer ester.
    nodes = {0: 'C', 1: 'O', 2: 'O', **{u: 'C' for u in range(3, 7)}}
    edges = [(0, 1, '='), (0, 2, '-'), (2, 3, '-'), (3, 4, '-'),
             (4, 5, '-'), (5, 6, '-'), (6, 0, '-')]
    assert not matches(chemistry.lactone_core("B"), complete_with_h(make_graph(nodes, edges)))


def test_rdkit_accepts_ordinary_covalent_substrates():
    Chem = pytest.importorskip('rdkit.Chem')
    for graph in (chemistry.monomer("B"), chemistry.chain(("B",) * 0), chemistry.chain(("B",) * 2)):
        mol = to_rdkit(graph)
        assert mol.GetNumAtoms() == len(node_rows(graph))
        assert all(atom.GetNumRadicalElectrons() == 0 for atom in mol.GetAtoms())
        assert all(bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                                         Chem.BondType.AROMATIC) for bond in mol.GetBonds())


def test_molfoundry_recognizes_substrates_as_molecules():
    from mfobs.core.gml import graph_gml
    from mfobs.core.engine import mf
    for graph in (chemistry.monomer("B"), chemistry.chain(("B",) * 0), chemistry.chain(("B",) * 2)):
        assert mf.Graph.fromGMLString(graph_gml(graph), add=False).isMolecule


def test_public_direct_application_of_ring_opening_without_composition():
    from mfobs.core.gml import graph_gml, rule_gml
    import molfoundry as mf
    rule = mf.Rule.fromGMLString(rule_gml(make_rule('ROP-net', *chemistry.propagation_sides("B"))), add=False)
    for n in (0, 1):
        hosts = [mf.Graph.fromGMLString(graph_gml(g), add=False)
                 for g in (chemistry.chain(("B",) * n), chemistry.monomer("B"))]
        expected = mf.Graph.fromGMLString(graph_gml(chemistry.chain(("B",) * (n+1))), add=False)
        dg = mf.DG(graphDatabase=hosts)
        with dg.build() as builder:
            outcomes = builder.apply(hosts, rule, mode=mf.ApplyMode.UseAny)
            assert len(outcomes) == 1
            assert len(outcomes[0].targets) == 1
            assert outcomes[0].targets[0].graph.isomorphism(expected)
