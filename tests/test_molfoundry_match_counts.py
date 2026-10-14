"""Minimal public-API count check, independent of mfobs model backends.

Run: python -m pytest -q tests/test_molfoundry_match_counts.py
"""
from mfobs.core.graph import node_rows, make_graph, edge_rows, node_labels
from mfobs.core.rules import make_rule
from itertools import combinations, permutations, product
from math import perm
import random

import molfoundry as mf
import pytest


HOST = '''node [ id 0 label "A" ] node [ id 1 label "A" ]
node [ id 2 label "A" ] node [ id 3 label "A" ]
edge [ source 0 target 1 label "-" ]
edge [ source 1 target 2 label "-" ]
edge [ source 2 target 3 label "-" ]'''


def query(size):
    nodes = ' '.join(f'node [ id {i} label "A" ]' for i in range(size))
    return mf.Rule.fromGMLString('rule [ context [ ' + nodes + ' ] ]', add=False)


@pytest.mark.parametrize('size', [1, 2, 3])
def test_injective_count_agrees_with_public_full_overlap_composition(size):
    from mfobs.core.graph import observable

    expected = len(list(permutations(range(4), size)))
    binding = mf.Rule.fromGMLString('rule [ right [ ' + HOST + ' ] ]', add=False)
    products = mf.RCEvaluator().eval(
        binding * mf.rcCommon(connected=False, includeEmpty=False) * query(size),
        onlyUnique=False,
    )
    # A full match leaves no unmatched input vertices in the composite.
    assert sum(not any(v.leftLabel is not None for v in p.vertices) for p in products) == expected
    pattern = make_graph({i: 'A' for i in range(size)})
    host = make_graph({i: 'A' for i in range(4)}, [(i, i+1, '-') for i in range(3)])
    assert observable(pattern, host) == expected


@pytest.mark.parametrize('size', [1, 2, 3])
def test_public_dg_count_for_a_single_host(size):
    host = mf.Graph.fromGMLString('graph [ ' + HOST + ' ]', add=False)
    rule = query(size)
    dg = mf.DG(graphDatabase=[host])
    with dg.build() as builder:
        edges = builder.apply([host], rule, mode=mf.ApplyMode.UseAny)
    # UseAny uses the listed host once, so every derivation has it as its one
    # source and the counts are embeddings into that host.
    assert all(edge.numSources == 1 for edge in edges)
    actual = sum(edge.numMatches(rule) for edge in edges)
    assert actual == len(list(permutations(range(4), size)))


def public_count(rule, host):
    """Counts on exactly one supplied host, not extra DG species copies."""
    from mfobs.core.gml import graph_gml
    graphs = [mf.Graph.fromGMLString(graph_gml(host), add=False)] if node_rows(host) else []
    dg = mf.DG(graphDatabase=graphs)
    with dg.build() as builder:
        edges = builder.apply(graphs, rule, mode=mf.ApplyMode.UseAny)
    return sum(edge.numMatches(rule) for edge in edges
               if edge.numSources == int(any(v.leftLabel is not None for v in rule.vertices)))


def pattern_count(pattern, host, no_edges=(), adjacency=()):

    from mfobs.core.gml import rule_gml
    rule = mf.Rule.fromGMLString(rule_gml(make_rule(
        'query', pattern, pattern, tuple(no_edges), tuple(adjacency))), add=False)
    return public_count(rule, host)


def exhaustive_count(pattern, host, no_edges=(), adjacency=()):
    """Small test oracle: try every injection, without any graph-search pruning."""
    import operator
    operators = {'=': operator.eq, '<': operator.lt, '<=': operator.le,
                 '>': operator.gt, '>=': operator.ge}
    pl, hl = node_labels(pattern), node_labels(host)
    he = {(u, v): label for u, v, label in edge_rows(host)}
    total = 0
    for images in permutations(hl, len(pl)):
        mapping = dict(zip(pl, images))
        if any(pl[u] != hl[v] for u, v in mapping.items()):
            continue
        if any(he.get(tuple(sorted((mapping[u], mapping[v])))) != label
               for u, v, label in edge_rows(pattern)):
            continue
        if any(tuple(sorted((mapping[u], mapping[v]))) in he for u, v in no_edges):
            continue
        permitted = True
        for guard in adjacency:
            image = mapping[guard.center]
            neighbours = sum(bond in guard.edge_labels and
                             hl[v if u == image else u] in guard.node_labels
                             for u, v, bond in edge_rows(host) if image in (u, v))
            permitted &= operators[guard.op](neighbours, guard.count)
        total += permitted
    return total


@pytest.mark.parametrize('host_mask', range(8))
def test_every_three_vertex_edge_set_and_noedge_guard(host_mask):


    pairs = list(combinations(range(3), 2))
    host = make_graph({i: 'A' for i in range(3)},
                      [(u, v, '-') for j, (u, v) in enumerate(pairs) if host_mask >> j & 1])
    for positive_mask, negative_mask in product(range(8), repeat=2):
        pattern = make_graph({i: 'A' for i in range(3)},
                             [(u, v, '-') for j, (u, v) in enumerate(pairs) if positive_mask >> j & 1])
        negative = tuple(pair for j, pair in enumerate(pairs) if negative_mask >> j & 1)
        assert pattern_count(pattern, host, negative) == exhaustive_count(pattern, host, negative)


def test_random_labelled_disconnected_patterns_and_full_host_neighbours():
    from mfobs.core.conditions import Adjacency


    rng = random.Random(20261002)
    for trial in range(250):
        hn, pn = rng.randrange(1, 7), rng.randrange(1, 5)
        host = make_graph({i: rng.choice('AAB') for i in range(hn)},
                          [(u, v, rng.choice(['-', '='])) for u, v in combinations(range(hn), 2)
                           if rng.random() < .4])
        pattern = make_graph({i: rng.choice('AAB') for i in range(pn)},
                             [(u, v, rng.choice(['-', '='])) for u, v in combinations(range(pn), 2)
                              if rng.random() < .3])
        negative = tuple((u, v) for u, v in combinations(range(pn), 2) if rng.random() < .15)
        adjacency = tuple(Adjacency(i, (rng.choice('AB'),), (rng.choice(['-', '=']),),
                                    rng.choice(['=', '<', '>', '<=', '>=']), rng.randrange(4))
                          for i in range(pn) if rng.random() < .25)
        assert pattern_count(pattern, host, negative, adjacency) == exhaustive_count(
            pattern, host, negative, adjacency), trial


@pytest.mark.parametrize('n', range(1, 7))
def test_falling_factorials_include_every_orientation_but_no_extra_host(n):

    host = make_graph({i: 'A' for i in range(n)}, [(i, i+1, '-') for i in range(n-1)])
    for k in range(1, 5):
        assert pattern_count(make_graph({i: 'A' for i in range(k)}), host) == perm(n, k)


def test_empty_hosts_and_reindexing():

    empty = make_graph({})
    host = make_graph({11: 'A', 20: 'A', 99: 'A'}, [(11, 20, '-')])
    pattern = make_graph({4: 'A', 9: 'A'})
    assert pattern_count(pattern, empty) == 0
    assert pattern_count(pattern, host) == 6


def test_empty_identity_is_the_constant_one_also_as_a_dg_derivation():
    from mfobs.core.graph import observable
    empty, host = make_graph({}), make_graph({0: 'A'})
    # The empty match uses no graph: DG.apply derives the identity on the empty
    # graph, a hyperedge without sources or targets, with one match.
    for state in (empty, host):
        assert pattern_count(empty, state) == 1
        assert observable(empty, state) == 1


def test_nested_condition_checks_unmatched_context():

    rule = mf.Rule.fromGMLString('''rule [ context [ node [ id 0 label "C" ] ]
      condition [ not [ exists [ node [ id 1 label "O" ]
        edge [ source 0 target 1 label "-" ]
        not [ exists [ node [ id 2 label "H" ] edge [ source 1 target 2 label "-" ] ] ]
      ] ] ] ]''', add=False)
    for attached in product([False, True], repeat=4):
        possible = [(0, 2, '-'), (1, 3, '-'), (2, 4, '-'), (3, 5, '-')]
        host = make_graph({0: 'C', 1: 'C', 2: 'O', 3: 'O', 4: 'H', 5: 'H'},
                          [edge for edge, present in zip(possible, attached) if present])
        # Each C is accepted unless it has an O neighbour without H.
        assert public_count(rule, host) == sum(not attached[i] or attached[i+2] for i in (0, 1))


def test_identical_reactants_count_ordered_copies_and_products():
    from collections import Counter
    from mfobs.core.gml import graph_gml

    rule = mf.Rule.fromGMLString('''rule [
      left [ node [ id 0 label "A" ] node [ id 1 label "A" ] ]
      right [ node [ id 0 label "B" ] node [ id 1 label "B" ] ] ]''', add=False)
    for n in range(2, 6):
        host = make_graph({i: 'A' for i in range(n)})
        assert public_count(rule, host) == n*(n-1)
    atom = mf.Graph.fromGMLString('graph [ node [ id 0 label "A" ] ]', add=False)
    dg = mf.DG(graphDatabase=[atom])
    with dg.build() as builder:
        edge, = builder.apply([atom, atom], rule, mode=mf.ApplyMode.UseAny)
    assert edge.numSources == 2 and edge.numMatches(rule) == 2
    # On A-A-A-A, the B sites have distances 1, 2 or 3; count both rule orientations.
    host = mf.Graph.fromGMLString(graph_gml(make_graph({i: 'A' for i in range(4)},
        [(0, 1, '-'), (1, 2, '-'), (2, 3, '-')])), add=False)
    dg = mf.DG(graphDatabase=[host])
    with dg.build() as builder:
        edges = builder.apply([host], rule, mode=mf.ApplyMode.UseAny)
    import networkx as nx
    by_distance = Counter()
    for edge in edges:
        if edge.numSources != 1:
            continue
        target, = edge.targets
        graph = target.graph
        nx_graph = nx.Graph((e.source.id, e.target.id) for e in graph.edges)
        b0, b1 = [v.id for v in graph.vertices if v.stringLabel == 'B']
        by_distance[nx.shortest_path_length(nx_graph, b0, b1)] += edge.numMatches(rule)
    assert by_distance == {1: 6, 2: 4, 3: 2}


def test_dangling_and_rule_specific_counts():

    from mfobs.core.gml import graph_gml
    delete = mf.Rule.fromGMLString('rule [ left [ node [ id 0 label "A" ] ] ]', add=False)
    relabel = mf.Rule.fromGMLString('''rule [ left [ node [ id 0 label "A" ] ]
        right [ node [ id 0 label "B" ] ] ]''', add=False)
    for n in range(2, 7):
        host = make_graph({i: 'A' for i in range(n)}, [(i, i+1, '-') for i in range(n-2)])
        assert public_count(delete, host) == (2 if n == 2 else 1)
        assert public_count(relabel, host) == n
    host = mf.Graph.fromGMLString(graph_gml(make_graph({0: 'A', 1: 'A'})), add=False)
    dg = mf.DG(graphDatabase=[host])
    with dg.build() as builder:
        edges = builder.apply([host], delete, mode=mf.ApplyMode.UseAny)
        for edge in edges:
            assert edge.numMatches(relabel) == 0
        again = builder.apply([host], delete, mode=mf.ApplyMode.UseAny)
    assert sum(e.numMatches(delete) for e in edges) == sum(e.numMatches(delete) for e in again) == 2
