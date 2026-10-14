"""Chemical validity and independent direct-application guard diagnostics."""
from mfobs.core.molecules import formula as molecular_formula
from mfobs.core.gml import bond_order_sums
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import edge_rows, graph_key, has_edge, node_labels, node_rows
from mfobs.core.rules import rule_graph
from collections import Counter
from math import comb

import pytest

from mfobs.core.conditions import matches, transitions

from mfobs.models import nonadjacent_chemistry as chemistry


def test_ordinary_covalent_hosts_and_atom_balance():
    expected = [({'C': 2, 'H': 4}, chemistry.ethylene()),
                ({'C': 4, 'H': 6}, chemistry.conjugated_diene()),
                ({'C': 6, 'H': 10}, chemistry.separated_diene())]
    for formula, graph in expected:
        assert molecular_formula(graph) == formula
        assert set(node_labels(graph).values()) <= {'C', 'H'}
        assert {label for _, _, label in edge_rows(graph)} <= {'-', '='}
        sums = bond_order_sums(graph)
        assert all(sums[u] == (4 if label == 'C' else 1) for u, label in node_rows(graph))
        assert to_rdkit(graph).GetNumAtoms() == graph.numVertices
    for rule in chemistry.rules().values():
        assert node_rows(rule_graph(rule, "left")) == node_rows(rule_graph(rule, "right"))
        assert bond_order_sums(rule_graph(rule, "left")) == bond_order_sums(rule_graph(rule, "right"))


@pytest.mark.parametrize('n', [0, 1, 2, 3, 4])
def test_eight_atom_maps_per_unordered_ethylene_pair(n):
    rule = chemistry.guarded_rule()
    host = chemistry.ethylene_mixture(n)
    targets = transitions(rule, host)
    assert len(matches(rule, host)) == 8*comb(n, 2)
    assert sum(targets.values()) == 8*comb(n, 2)
    assert chemistry.PER_MATCH * sum(targets.values()) == comb(n, 2)
    for target in targets:
        assert node_rows(target) == node_rows(host)
        assert bond_order_sums(target) == bond_order_sums(host)
        assert to_rdkit(target).GetNumAtoms() == host.numVertices


def test_context_guard_is_genuinely_active_with_valid_bicyclobutane_ablation():
    host = chemistry.conjugated_diene()
    assert has_edge(host, 0, 2, '-')
    assert not transitions(chemistry.guarded_rule(), host)
    ablated = transitions(chemistry.ablated_rule(), host)
    assert sum(ablated.values()) == 4
    for target in ablated:
        assert node_rows(target) == node_rows(host)
        assert bond_order_sums(target) == bond_order_sums(host)
        heavy = [(u, v, label) for u, v, label in edge_rows(target)
                 if node_labels(target)[u] == node_labels(target)[v] == 'C']
        assert len(heavy) == 5
        assert all(label == '-' for _, _, label in heavy)
        degrees = Counter(vertex for u, v, _ in heavy for vertex in (u, v))
        assert sorted(degrees.values()) == [2, 2, 3, 3]
        assert molecular_formula(target) == {'C': 4, 'H': 6}
        assert to_rdkit(target).GetNumAtoms() == 10
        # The pre-existing context bond is kept, not created by the rule.
        assert has_edge(target, 0, 2, '-')


def test_separated_diene_is_not_excluded_by_context_guards():
    host = chemistry.separated_diene()
    guarded = transitions(chemistry.guarded_rule(), host)
    ablated = transitions(chemistry.ablated_rule(), host)
    assert sum(guarded.values()) == 8
    assert {graph_key(g): n for g, n in guarded.items()} == {
        graph_key(g): n for g, n in ablated.items()}
    for target in guarded:
        assert to_rdkit(target).GetNumAtoms() == host.numVertices
        assert bond_order_sums(target) == bond_order_sums(host)


@pytest.mark.parametrize('n,p', [(0, 2), (1, 2), (2, 1), (3, 2), (4, 2)])
def test_depletion_sublanguage_preserved_with_inert_products(n, p):
    from rdkit import Chem
    def key(graph):
        return Chem.MolToSmiles(Chem.RemoveHs(to_rdkit(graph)))
    host = chemistry.disjoint_union([chemistry.ethylene()]*n + [chemistry.cyclobutane()]*p)
    targets = transitions(chemistry.guarded_rule(), host)
    assert sum(targets.values()) == 8*comb(n, 2)
    if n >= 2:
        expected = chemistry.disjoint_union([chemistry.ethylene()]*(n-2) +
                                           [chemistry.cyclobutane()]*(p+1))
        assert all(key(target) == key(expected) for target in targets)
    assert not transitions(chemistry.guarded_rule(), chemistry.cyclobutane())


def test_public_direct_applications_agree_with_independent_successors():
    from mfobs.core.gml import component_graphs, graph_gml, rule_gml
    import molfoundry as mf
    def keys(graph):
        return tuple(sorted(mf.Graph.fromGMLString(graph_gml(part), add=False).canonKey
                            for part in component_graphs(graph)))
    hosts = [chemistry.ethylene_mixture(1), chemistry.ethylene_mixture(2),
             chemistry.conjugated_diene(), chemistry.separated_diene()]
    checked = 0
    observed = []
    for host in hosts:
        parts = [mf.Graph.fromGMLString(graph_gml(part), add=False)
                 for part in component_graphs(host)]
        initial = Counter(part.canonKey for part in parts)
        species = list({part.canonKey: part for part in parts}.values())
        for name, rule in chemistry.rules().items():
            engine_rule = mf.Rule.fromGMLString(rule_gml(rule), add=False)
            dg = mf.DG(graphDatabase=species)
            actual = set()
            with dg.build() as builder:
                for edge in builder.apply(species, engine_rule, mode=mf.ApplyMode.UseAnyMultiplicity):
                    # DG.apply explores species derivations and can reuse a
                    # supplied species twice. The finite-host oracle instead
                    # has a finite molecular inventory: enforce it explicitly.
                    consumed = Counter(v.graph.canonKey for v in edge.sources)
                    if any(amount > initial[key] for key, amount in consumed.items()):
                        continue
                    remainder = initial.copy()
                    remainder.subtract(consumed)
                    remainder.update(v.graph.canonKey for v in edge.targets)
                    actual.add(tuple(sorted(remainder.elements())))
            observed.append((name, actual, rule, host))
            checked += 1
    assert checked == 8
    for label, actual, rule, host in observed:
        assert actual == {keys(target) for target in transitions(rule, host)}, label
