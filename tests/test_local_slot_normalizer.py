"""Actual symbolic normalizer, checked independently without finite-basis input."""
from mfobs.core.graph import make_graph, edge_rows, node_labels
from mfobs.core.rules import make_rule, is_identity, rule_graph
from itertools import combinations, product

import pytest

from mfobs.core.conditions import Adjacency, parse_observable, engine_no_edges
from mfobs.core.engine import mf
from mfobs.core.gml import rule_gml

from mfobs.invariants.local_slots import CONTRACT, Normalizer, meaning
from test_local_slot_preservation import (
    count, current_rules, identity, normalization_vectors, slot_host,
)


def normalized_count(reduction, host):
    assert reduction.status == 'exact'
    pattern, adj = parse_observable(reduction.rule)
    literal = identity(make_graph(pattern.labels, pattern.edges), pattern.no_edges, adj)
    return count(literal, host)


def small_hosts():
    yield slot_host(())
    for n in range(1, 4):
        topologies = [(), tuple(combinations(range(n), 2))]
        if n == 3:
            topologies.append(((0, 1), (1, 2)))
        for edges in set(topologies):
            for bits in product((0, 1), repeat=n):
                yield slot_host(bits, edges)
    for bits in ((0, 0, 0, 0), (1, 0, 1, 0), (1, 1, 1, 1)):
        yield slot_host(bits, ((0, 1), (0, 2), (0, 3)))
    # The encoded local CONTRACT permits unowned isolated P/L vertices
    # (at most one owner), unlike the stricter complete-host sample generator.
    # Complete anchored patterns must remain insensitive to those extra ports.
    host = slot_host((1, 0), ((0, 1),))
    yield make_graph({**node_labels(host), 90: 'P', 91: 'L'}, edge_rows(host))


@pytest.mark.parametrize('name', [
    'shared_P', 'shared_L', 'multiple_P', 'multiple_L',
    'cross_owner_bound', 'endpoint_two_partners',
])
def test_local_positive_contradictions_are_zero_before_completion(name):
    vector = next(item for item in normalization_vectors() if item['name'] == name)
    reduction = Normalizer().normalize(vector['rule'])
    assert reduction.status == 'zero'
    assert reduction.rule is None and reduction.key is None
    assert all(count(vector['rule'], host) == 0 for host in small_hosts())


@pytest.mark.parametrize('graph', [
    make_graph({0: 'P', 1: 'P'}, [(0, 1, 'bound')]),
    make_graph({0: 'C', 1: 'C'}, [(0, 1, 'cp')]),
    make_graph({0: 'P', 1: 'L'}, [(0, 1, 'next')]),
    make_graph({0: 'C', 1: 'P'}, [(0, 1, 'not-in-signature')]),
    make_graph({0: 'X'}),
])
def test_forbidden_signature_is_zero(graph):
    assert Normalizer().normalize_graph(graph).status == 'zero'


@pytest.mark.parametrize('graph,no_edges', [
    (make_graph({0: 'C'}), ()),
    (make_graph({0: 'P'}), ()),
    (make_graph({0: 'L'}), ()),
    (make_graph({0: 'C', 1: 'P'}, [(0, 1, 'cp')]), ()),
    (make_graph({0: 'P', 1: 'L'}, [(0, 1, 'bound')]), ()),
    (make_graph({0: 'P', 1: 'L'}), ((0, 1),)),
    (make_graph({0: 'C', 1: 'P', 2: 'L'}, [(0, 1, 'cp'), (1, 2, 'bound')]), ((0, 2),)),
])
def test_partial_ownership_is_unknown_not_assumed_completed(graph, no_edges):
    reduction = Normalizer().normalize(identity(graph, no_edges))
    assert reduction.status == 'unknown'
    assert 'partial ownership' in reduction.reason


def test_unknown_partial_P_L_NAC_includes_cross_owner_matches():
    rule = identity(make_graph({0: 'P', 1: 'L'}), ((0, 1),))
    host = slot_host((1, 0))
    assert count(rule, host) == 3
    assert Normalizer().normalize(rule).status == 'unknown'
    assert count(identity(slot_host((0,)), ((1, 2),)), host) == 1


def test_positive_and_negative_edge_contradiction_is_zero():
    for graph, no_edges in [(slot_host((1,)), ((1, 2),)),
                            (slot_host((0,)), ((0, 1),))]:
        assert Normalizer().normalize_graph(graph, no_edges).status == 'zero'


def test_same_owner_NAC_retained_cross_owner_and_wrong_type_NACs_removed():
    graph = slot_host((0, 0))
    same_owner = ((1, 2),)
    implied_absent = ((1, 5), (2, 4), (1, 4), (2, 5), (0, 4), (0, 5))
    rule = identity(graph, same_owner+implied_absent)
    reduction = Normalizer().normalize_graph(graph, engine_no_edges(rule))
    assert reduction.status == 'exact'
    assert reduction.quotient['occupancies'] == [0, None]
    assert set(reduction.removed_no_edges) == set(implied_absent)
    _, adjacency = parse_observable(reduction.rule)
    assert not adjacency
    assert len(reduction.rule.matchConstraintsNoEdge) == 1
    for host in small_hosts():
        assert normalized_count(reduction, host) == count(rule, host)


def test_center_NoEdge_is_not_removed_for_arbitrary_next_topology():
    rule = identity(slot_host((0, 0)), ((0, 3),))
    reduction = Normalizer().normalize(rule)
    assert reduction.status == 'exact'
    assert reduction.quotient['no_next_edges'] == [[0, 1]]
    assert not reduction.removed_no_edges
    assert normalized_count(reduction, slot_host((0, 0))) == 2
    assert normalized_count(reduction, slot_host((0, 0), ((0, 1),))) == 0


def test_wildcard_free_and_occupied_are_three_different_observables():
    normalizer = Normalizer()
    wildcard = normalizer.normalize(identity(slot_host((0,))))
    free = normalizer.normalize(identity(slot_host((0,)), ((1, 2),)))
    occupied = normalizer.normalize(identity(slot_host((1,))))
    assert wildcard.quotient['occupancies'] == [None]
    assert free.quotient['occupancies'] == [0]
    assert occupied.quotient['occupancies'] == [1]
    assert len({reduction.key for reduction in (wildcard, free, occupied)}) == 3
    assert meaning(wildcard.quotient) == '* (unspecified)'
    assert meaning(free.quotient) == '0 (free)'
    host = slot_host((0, 1, 1))
    assert [normalized_count(item, host) for item in (wildcard, free, occupied)] == [3, 1, 2]


def test_empty_pattern_is_exact_constant_including_empty_host():
    reduction = Normalizer().normalize(identity(make_graph({})))
    assert reduction.status == 'exact'
    assert reduction.quotient['site_count'] == 0
    assert meaning(reduction.quotient) == 'empty-pattern constant 1'
    assert all(normalized_count(reduction, host) == 1 for host in small_hosts())


def test_counts_preserved_without_automorphism_division_on_anchored_fragment():
    normalizer, checked = Normalizer(), 0
    hosts = list(small_hosts())
    for n in range(4):
        possible = tuple(combinations(range(n), 2))
        topologies = [(), possible]
        if n == 3:
            topologies.append(((0, 1), (1, 2)))
        for occupancy in product((None, 0, 1), repeat=n):
            for edges in set(topologies):
                graph = slot_host(tuple(int(value == 1) for value in occupancy), edges)
                no_edges = [(3*i+1, 3*i+2) for i, value in enumerate(occupancy) if value == 0]
                no_edges += [(3*i+1, 3*j+2) for i in range(n) for j in range(n) if i != j]
                original = identity(graph, no_edges)
                reduction = normalizer.normalize(original)
                assert reduction.status == 'exact'
                assert reduction.record()['multiplicity_factor'] == 1
                for host in hosts:
                    assert normalized_count(reduction, host) == count(original, host)
                    checked += 1
    assert checked > 3000
    symmetric = identity(slot_host((0, 0), ((0, 1),)))
    normalized = normalizer.normalize(symmetric)
    assert normalized_count(normalized, slot_host((0, 0), ((0, 1),))) == 2


def test_diagonal_delete_recreate_bound_span_flattens_to_strict_identity():
    jump = mf.Rule.fromGMLString('''rule [ ruleID "diagonal-bound"
        context [ node [ id 0 label "C" ] node [ id 1 label "P" ] node [ id 2 label "L" ]
          edge [ source 0 target 1 label "cp" ] edge [ source 0 target 2 label "cl" ] ]
        left [ edge [ source 1 target 2 label "bound" ] ]
        right [ edge [ source 1 target 2 label "bound" ] ] ]''', add=False)
    reduction = Normalizer().normalize(jump)
    assert reduction.status == 'exact'
    assert is_identity(reduction.rule)
    for host in small_hosts():
        assert normalized_count(reduction, host) == count(identity(slot_host((1,))), host)


def test_non_diagonal_reactions_are_unknown_to_observable_normalizer():
    normalizer = Normalizer()
    for rule in current_rules().values():
        assert normalizer.normalize(rule).status == 'unknown'


def test_normalizer_cache_is_canonical_and_does_not_rescale_counts():
    rule = identity(slot_host((0,)), ((1, 2),))
    normalizer = Normalizer()
    first = normalizer.normalize(rule)
    assert normalizer.normalize(rule) is first
    assert normalizer.stats['exact'] == 1 and len(normalizer.cache) == 1
    assert first.key == first.rule.canonKey
    assert first.record()['multiplicity_factor'] == 1


def test_unknown_conditions_and_malformed_NoEdges_fail_closed():
    normalizer = Normalizer()
    slot = slot_host((0,))
    for no_edges in (((1, 90),), ((1, 1),), ((1,),)):
        assert normalizer.normalize_graph(slot, no_edges).status == 'unknown'
    adjacency = (Adjacency(0, ('C',), ('next',), '<', 4),)
    assert normalizer.normalize(identity(slot, adjacency=adjacency)).status == 'unknown'


def test_degree_bound_is_optional_and_must_be_supplied_explicitly():
    graph = slot_host((0, 0, 0, 0), ((0, 1), (0, 2), (0, 3)))
    assert Normalizer().normalize_graph(graph).status == 'exact'
    bounded = Normalizer({**CONTRACT, 'next_max_degree': 2})
    assert bounded.normalize_graph(graph).status == 'zero'
    for value in (-1, True, 2.5):
        with pytest.raises(ValueError):
            Normalizer({**CONTRACT, 'next_max_degree': value})


def test_contract_changes_outside_reviewed_fragment_are_rejected():
    with pytest.raises(ValueError):
        Normalizer({**CONTRACT, 'bound_endpoints_have_the_same_C': False})
    with pytest.raises(ValueError):
        Normalizer({**CONTRACT, 'sites': 8})


def test_current_three_rules_have_syntactic_preservation_certificates():
    normalizer = Normalizer()
    for name, rule in current_rules().items():
        certificate = normalizer.certify_rule_preservation(rule)
        assert certificate['passed']
        assert certificate['vertices_preserved'] and certificate['ownership_and_next_preserved']
        assert len(certificate['added_bound_edges']) == int(name != 'unbind')
        assert len(certificate['deleted_bound_edges']) == int(name == 'unbind')
        assert certificate['method'] == 'syntactic sufficient-condition check, not sampled-state testing'


@pytest.mark.parametrize('name', ['basal_bind', 'cooperative_bind'])
def test_molfoundry_injects_created_edge_NoEdge_omitted_from_literal_GML(name):
    source = current_rules()[name]
    rule = make_rule(source.name, rule_graph(source), rule_graph(source, "right"))
    gml = rule_gml(rule)
    assert 'constrainNoEdge' not in gml
    imported = mf.Rule.fromGMLString(gml, add=False)
    # The public importer adds the simple-graph edge-creation safeguard.
    # Literal omission therefore is NOT an ablation of the effective rule.
    assert imported.matchConstraintsNoEdge
    assert all(c.implied for c in imported.matchConstraintsNoEdge)
    assert Normalizer().certify_rule_preservation(imported)['passed']


@pytest.mark.parametrize('name', ['basal_bind', 'cooperative_bind'])
def test_checker_rejects_a_genuinely_unguarded_effective_rule_record(name, monkeypatch):
    import mfobs.invariants.local_slots as local_slots
    imported = current_rules()[name]
    # A negative unit test of the sufficient-condition checker. No public
    # engine object is changed; the ordinary importer repairs this omission.
    monkeypatch.setattr(local_slots, 'engine_no_edges', lambda _: ())
    with pytest.raises(ValueError, match='explicit absence guard'):
        Normalizer().certify_rule_preservation(imported)


def test_cross_owner_binding_cannot_be_certified_even_with_target_NoEdge():
    left = slot_host((0, 0))
    right = make_graph(node_labels(left), [*edge_rows(left), (1, 5, 'bound')])
    rule = make_rule('cross-owner-bind', left, right, ((1, 5),))
    with pytest.raises(ValueError, match='different known owners'):
        Normalizer().certify_rule_preservation(rule)


def test_rule_topology_or_vertex_mutations_are_not_covered():
    left = slot_host((0, 0), ((0, 1),))
    right = make_graph(node_labels(left), [edge for edge in edge_rows(left) if edge[2] != 'next'])
    with pytest.raises(ValueError, match='topology changes'):
        Normalizer().certify_rule_preservation(make_rule('delete-next', left, right))
    right = slot_host((0,))
    with pytest.raises(ValueError, match='vertices to be preserved'):
        Normalizer().certify_rule_preservation(make_rule('delete-slot', left, right))
