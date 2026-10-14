"""Independent local-slot invariant tests; finite hosts are VERIFICATION ONLY.

Run from the molfoundry-observables project with PYTHONPATH=src. These tests
import neither the proposed local_slots normalizer nor a synthesis adapter.
``normalization_vectors`` exports positive/negative test patterns for the
eventual Counter.normalize API. No finite host table is supplied to synthesis.

The main host corpus has arbitrary finite SIMPLE undirected ``next`` topology.
Every C owns exactly one P by cp and one L by cl; every P and L here has
exactly one C owner. This fully-owned subset is STRONGER than production
CONTRACT's at-most-one-owner premise, which also permits inert orphan ports.
Those extra ports are checked separately below and by the normalizer tests.
A bound edge is allowed only between its owner's P/L pair. Endpoint bound
capacity <=1 follows from these ownership assumptions, but is also checked
explicitly to expose weakened-assumption counterexamples.
"""
from mfobs.core.conditions import engine_conditions, engine_no_edges, matches, transitions
from mfobs.core.graph import node_rows, graph_key, make_graph, edge_rows, node_labels
from mfobs.core.rules import make_rule, rule_graph
from collections import Counter
from itertools import combinations, product
import os
from pathlib import Path
import random

import pytest


from mfobs.core.gml import read_rule_record


PROJECT = Path(os.environ.get('MF_LOCAL_SLOT_PROJECT', Path(__file__).resolve().parents[1]))
ALLOWED = {'cp': frozenset(('C', 'P')), 'cl': frozenset(('C', 'L')),
           'bound': frozenset(('P', 'L')), 'next': frozenset(('C',))}


def slot_host(bits, next_edges=()):
    """Stable IDs (C,P,L)=(3*i,3*i+1,3*i+2); no topology restrictions."""
    bits = tuple(bits)
    if any(type(bit) is not int or bit not in (0, 1) for bit in bits):
        raise ValueError('bits must be integers 0 or 1')
    nodes, edges = {}, []
    for i, bit in enumerate(bits):
        c, p, l = 3*i, 3*i+1, 3*i+2
        nodes.update({c: 'C', p: 'P', l: 'L'})
        edges += [(c, p, 'cp'), (c, l, 'cl')]
        if bit:
            edges.append((p, l, 'bound'))
    for i, j in next_edges:
        if not (0 <= i < len(bits) and 0 <= j < len(bits)) or i == j:
            raise ValueError('next topology must be simple and use existing slots')
        edges.append((3*i, 3*j, 'next'))
    return make_graph(nodes, edges)


def identity(graph, no_edges=(), adjacency=(), name='test-observable'):
    return make_rule(name, graph, graph, tuple(no_edges), tuple(adjacency))


def count(rule, host):
    return len(matches(rule_graph(rule), host, engine_no_edges(rule), engine_conditions(rule)))


def ownership_violations(graph):
    """Independent full-host validator, not an embedding-pattern validator."""
    labels, failures = node_labels(graph), []
    if set(labels.values())-{'C', 'P', 'L'}:
        failures.append('unsupported vertex type')
    for u, v, edge in edge_rows(graph):
        if edge not in ALLOWED or frozenset((labels[u], labels[v])) != ALLOWED[edge]:
            failures.append('bad edge typing')
    def around(u, edge, kind):
        return [b if a == u else a for a, b, label in edge_rows(graph)
                if u in (a, b) and label == edge and labels[b if a == u else a] == kind]
    owners = {}
    for u, kind in node_rows(graph):
        if kind == 'C':
            for edge, port in (('cp', 'P'), ('cl', 'L')):
                if len(around(u, edge, port)) != 1:
                    failures.append('C must have exactly one '+port)
        if kind in ('P', 'L'):
            candidates = around(u, 'cp' if kind == 'P' else 'cl', 'C')
            if len(candidates) != 1:
                failures.append(kind+' must have exactly one owner')
            else:
                owners[u] = candidates[0]
            if sum(u in (a, b) and edge == 'bound' for a, b, edge in edge_rows(graph)) > 1:
                failures.append('bound endpoint capacity exceeded')
    for u, v, edge in edge_rows(graph):
        if edge == 'bound' and (u not in owners or v not in owners or owners[u] != owners[v]):
            failures.append('cross-owner or unowned bound')
    return failures


def current_rules():
    return {name: read_rule_record((PROJECT/'examples/ring8/rules'/f'{name}.gml').read_text(encoding='utf-8'))
            for name in ('basal_bind', 'cooperative_bind', 'unbind')}


def verification_hosts():
    """All topologies/occupancies through four sites, plus seeded larger tests."""
    for n in range(5):
        possible = tuple(combinations(range(n), 2))
        for topology_bits in product((0, 1), repeat=len(possible)):
            edges = tuple(edge for edge, bit in zip(possible, topology_bits) if bit)
            for bits in product((0, 1), repeat=n):
                yield bits, edges, slot_host(bits, edges)
    rng = random.Random(20260909)
    for n in range(5, 9):
        possible = tuple(combinations(range(n), 2))
        special = [(), tuple((i, i+1) for i in range(n-1)),
                   tuple((0, i) for i in range(1, n)), possible,
                   tuple((i, (i+1) % n) for i in range(n))]
        for _ in range(15):
            special.append(tuple(edge for edge in possible if rng.random() < .35))
        for edges in special:
            bits = tuple(rng.randrange(2) for _ in range(n))
            yield bits, edges, slot_host(bits, edges)


def verify_preservation():
    rules = current_rules()
    summary = {'hosts': 0, 'source_rule_checks': 0, 'raw_matches': 0,
               'distinct_labelled_successors': 0, 'maximum_next_degree': 0}
    for bits, edges, host in verification_hosts():
        assert not ownership_violations(host)
        degree = Counter(u for edge in edges for u in edge)
        summary['maximum_next_degree'] = max(summary['maximum_next_degree'], max(degree.values(), default=0))
        summary['hosts'] += 1
        for name, rule in rules.items():
            actual = transitions(rule, host)
            expected = Counter()
            for i, bit in enumerate(bits):
                multiplicity = (1 if bit == 0 else 0) if name == 'basal_bind' else (
                    (1 if bit else 0) if name == 'unbind' else
                    (sum(bits[b if a == i else a] for a, b in edges if i in (a, b)) if bit == 0 else 0))
                if multiplicity:
                    changed = list(bits)
                    changed[i] = 1-bit
                    expected[slot_host(changed, edges)] += multiplicity
            assert Counter({graph_key(g): n for g, n in actual.items()}) == Counter({graph_key(g): n for g, n in expected.items()}), (name, bits, edges)
            assert sum(actual.values()) == count(rule, host)
            summary['source_rule_checks'] += 1
            summary['raw_matches'] += sum(actual.values())
            summary['distinct_labelled_successors'] += len(actual)
            for target in actual:
                assert node_rows(target) == node_rows(host)
                assert not ownership_violations(target)
                assert {(u, v, edge) for u, v, edge in edge_rows(target) if edge != 'bound'} == {
                    (u, v, edge) for u, v, edge in edge_rows(host) if edge != 'bound'}
    return summary


def bad_hosts():
    base = slot_host((0, 0))
    nodes, edges = node_labels(base), list(edge_rows(base))
    return {
        'shared_P_owner': make_graph(nodes, edges+[(3, 1, 'cp')]),
        'shared_L_owner': make_graph(nodes, edges+[(3, 2, 'cl')]),
        'one_C_multiple_P': make_graph({**nodes, 6: 'P'}, edges+[(0, 6, 'cp')]),
        'one_C_multiple_L': make_graph({**nodes, 6: 'L'}, edges+[(0, 6, 'cl')]),
        'cross_owner_bound': make_graph(nodes, edges+[(1, 5, 'bound')]),
        'two_bound_partners': make_graph(nodes, edges+[(1, 2, 'bound'), (1, 5, 'bound')]),
        'swapped_partners_capacity_one': make_graph(nodes, edges+[(1, 5, 'bound'), (4, 2, 'bound')]),
        'missing_P': make_graph({u: kind for u, kind in nodes.items() if u != 1},
                                     [edge for edge in edges if 1 not in edge[:2]]),
        'bad_bound_typing': make_graph(nodes, edges+[(1, 4, 'bound')]),
    }


def normalization_vectors():
    """Counter.normalize test vectors without assuming its return schema.

    preferred_status=zero means zero on ALL admitted hosts, not just this corpus.
    Unknown is always a safe fallback. Exact cases include the original rule;
    a normalizer need not eagerly complete a partial ownership graph.
    """
    graph = make_graph
    vectors = [
        ('shared_P', graph({0:'C',1:'P',2:'C'}, [(0,1,'cp'),(2,1,'cp')]), (), 'zero'),
        ('shared_L', graph({0:'C',1:'L',2:'C'}, [(0,1,'cl'),(2,1,'cl')]), (), 'zero'),
        ('multiple_P', graph({0:'C',1:'P',2:'P'}, [(0,1,'cp'),(0,2,'cp')]), (), 'zero'),
        ('multiple_L', graph({0:'C',1:'L',2:'L'}, [(0,1,'cl'),(0,2,'cl')]), (), 'zero'),
        ('cross_owner_bound', graph({0:'C',1:'P',2:'C',3:'L'},
                                   [(0,1,'cp'),(2,3,'cl'),(1,3,'bound')]), (), 'zero'),
        ('endpoint_two_partners', graph({0:'P',1:'L',2:'L'}, [(0,1,'bound'),(0,2,'bound')]), (), 'zero'),
        ('hidden_ownership_edge_contradicts_NoEdge', graph({0:'C',1:'P',2:'L'},
                                   [(0,1,'cp'),(1,2,'bound')]), ((0,2),), 'zero'),
        ('complete_free_slot', slot_host((0,)), ((1,2),), 'exact'),
        ('complete_occupancy_unknown', slot_host((0,)), (), 'exact'),
        ('bare_C', graph({0:'C'}), (), 'exact'),
        ('bare_P', graph({0:'P'}), (), 'exact'),
        ('bare_L', graph({0:'L'}), (), 'exact'),
        ('bare_bound_pair', graph({0:'P',1:'L'}, [(0,1,'bound')]), (), 'exact'),
        ('partial_unbound_pair_ambiguous_owners', graph({0:'P',1:'L'}), ((0,1),), 'exact'),
        ('partial_P_L_ambiguous_owners', graph({0:'P',1:'L'}), (), 'exact'),
        ('P_P_NoEdge_redundant_by_typing', graph({0:'P',1:'P'}), ((0,1),), 'exact'),
        ('C_C_NoEdge_not_redundant', graph({0:'C',1:'C'}), ((0,1),), 'exact'),
        ('degree_three_next_star_is_valid', graph({u:'C' for u in range(4)},
                                   [(0,1,'next'),(0,2,'next'),(0,3,'next')]), (), 'exact'),
    ]
    return [{'name': name, 'rule': identity(g, ne, name=name), 'preferred_status': status}
            for name, g, ne, status in vectors]


def test_arbitrary_small_topology_rules_preserve_full_contract():
    result = verify_preservation()
    assert result['hosts'] > 1000
    assert result['maximum_next_degree'] == 7  # complete eight-site graph; NOT a ring assumption


@pytest.mark.parametrize('name', tuple(bad_hosts()))
def test_invalid_full_hosts_are_rejected(name):
    assert ownership_violations(bad_hosts()[name])


def test_capacity_one_alone_does_not_make_current_NoEdge_sufficient():
    host = bad_hosts()['swapped_partners_capacity_one']
    assert set(ownership_violations(host)) == {'cross-owner or unowned bound'}
    basal = current_rules()['basal_bind']
    products = transitions(basal, host)
    assert sum(products.values()) == 2
    assert all('bound endpoint capacity exceeded' in ownership_violations(target) for target in products)


def test_capacity_is_a_consequence_of_unique_own_partner_not_a_new_guard():
    for bits, edges, host in verification_hosts():
        for u, kind in node_rows(host):
            if kind in ('P', 'L'):
                bound_degree = sum(u in (a, b) and edge == 'bound' for a, b, edge in edge_rows(host))
                assert bound_degree in (0, 1)
                i = u//3
                assert bound_degree == bits[i]


def test_inert_orphan_ports_allowed_by_at_most_contract_do_not_change_rules():
    """Extend the fully-owned corpus toward the weaker production contract."""
    for bits in ((0, 0), (0, 1), (1, 0), (1, 1)):
        base = slot_host(bits, ((0, 1),))
        extended = make_graph({**node_labels(base), 90: 'P', 91: 'L'}, edge_rows(base))
        for rule in current_rules().values():
            expected = Counter({make_graph({**node_labels(target), 90: 'P', 91: 'L'}, edge_rows(target)): multiplicity
                                for target, multiplicity in transitions(rule, base).items()})
            assert Counter({graph_key(g): n for g, n in transitions(rule, extended).items()}) == Counter({graph_key(g): n for g, n in expected.items()})


def test_zero_vectors_have_no_admitted_matches():
    zero = [vector for vector in normalization_vectors() if vector['preferred_status'] == 'zero']
    for bits in product((0, 1), repeat=4):
        host = slot_host(bits, tuple(combinations(range(4), 2)))
        assert all(count(vector['rule'], host) == 0 for vector in zero)


def test_unique_owner_completions_preserve_counts_without_automorphism_division():
    partials = [make_graph({0:kind}) for kind in ('C','P','L')]
    partials += [make_graph({0:'C',1:'P'}, [(0,1,'cp')])]
    bound_pair = identity(make_graph({0:'P',1:'L'}, [(0,1,'bound')]))
    full = identity(slot_host((0,)))
    occupied = identity(slot_host((1,)))
    for n in range(5):
        for bits in product((0,1), repeat=n):
            host = slot_host(bits, tuple(combinations(range(n),2)))
            assert count(full, host) == n
            assert all(count(identity(pattern),host) == n for pattern in partials)
            assert count(bound_pair,host) == count(occupied,host) == sum(bits)
    two_unknown_slots = identity(slot_host((0,0), ((0,1),)))
    assert count(two_unknown_slots,slot_host((0,0),((0,1),))) == 2


def test_partial_P_L_pair_must_not_be_forced_to_share_owner():
    bare = make_graph({0:'P',1:'L'})
    unconstrained, negative = identity(bare), identity(bare,((0,1),))
    same_owner_free = identity(slot_host((0,)),((1,2),))
    for n in range(5):
        for bits in product((0,1), repeat=n):
            host = slot_host(bits)
            assert count(unconstrained,host) == n*n
            assert count(negative,host) == n*n-sum(bits)
            assert count(same_owner_free,host) == n-sum(bits)
    host = slot_host((1,0))
    assert count(negative,host) == 3 and count(same_owner_free,host) == 1


def test_next_absence_is_genuine_and_next_degree_is_not_two():
    nonadjacent = identity(make_graph({0:'C',1:'C'}),((0,1),))
    assert count(nonadjacent,slot_host((0,0))) == 2
    assert count(nonadjacent,slot_host((0,0),((0,1),))) == 0
    star = next(vector['rule'] for vector in normalization_vectors()
                if vector['name']=='degree_three_next_star_is_valid')
    assert count(star,slot_host((0,0,0,0),((0,1),(0,2),(0,3)))) == 6


if __name__ == '__main__':
    import json
    print(json.dumps(verify_preservation(), indent=2))
