"""Independent regression tests for joint injective counts on one chain."""
from mfobs.core.graph import make_graph
from math import comb, factorial
import random
import pytest

from mfobs.core.conditions import Adjacency, matches
from mfobs.invariants.chain import ChainLanguage, ChainCounter, CountingUnknown


def path_language():
    g = {'nodes': [[0, 'C']], 'edges': []}
    return ChainLanguage({'kind': 'one-port-decorated-path', 'join_bond': '-', 'zero': g,
        'left': {**g, 'incoming': None, 'outgoing': 0},
        'repeat': {**g, 'incoming': 0, 'outgoing': 0},
        'right': {**g, 'incoming': 0, 'outgoing': None}})


def test_arbitrarily_high_free_vertex_factorials():
    counter = ChainCounter(path_language())
    for k in [2, 4, 10, 20]:
        p = make_graph({i: 'C' for i in range(k)})
        result = counter.count(p)
        assert result.witness['method'] == 'free-vertex-falling-factorials'
        for n in [0, 1, 5, 20, 1000]:
            want = factorial(k) * comb(n + 1, k) if n + 1 >= k else 0
            assert result.sequence.value(n) == want


def test_marked_edge_and_free_vertices():
    counter = ChainCounter(path_language())
    for k in [1, 4, 10]:
        p = make_graph({i: 'C' for i in range(k + 2)}, [(0, 1, '-')])
        result = counter.count(p)
        for n in [0, 1, 5, 20, 1000]:
            want = 2*n*factorial(k)*comb(n - 1, k) if n - 1 >= k else 0
            assert result.sequence.value(n) == want


def test_two_disjoint_edges_kernel_correction():
    counter = ChainCounter(path_language())
    p = make_graph({i: 'C' for i in range(4)}, [(0, 1, '-'), (2, 3, '-')])
    result = counter.count(p)
    assert result.witness['method'] == 'component-kernel-partition'
    for n in [0, 1, 2, 3, 7, 1000]:
        assert result.sequence.value(n) == (4*(n-1)*(n-2) if n >= 2 else 0)


def test_three_isolated_nonadjacent_vertices():
    counter = ChainCounter(path_language())
    p = make_graph({i: 'C' for i in range(3)})
    result = counter.count(p, ((0, 1), (1, 2), (0, 2)))
    assert result.witness['method'] == 'cross-noedge-complement'
    for n in [0, 1, 2, 3, 4, 6, 10, 1000]:
        assert result.sequence.value(n) == (6*comb(n-1, 3) if n >= 4 else 0)


def test_conditioned_isolated_vertices_and_short_chain_exception():
    counter = ChainCounter(path_language())
    p = make_graph({0: 'C', 1: 'C'})
    conditions = tuple(Adjacency(i, ('C',), ('-',), '=', 1) for i in range(2))
    ordinary = counter.count(p, (), conditions)
    no_edge = counter.count(p, ((0, 1),), conditions)
    for n in [0, 1, 2, 3, 30, 1000]:
        assert ordinary.sequence.value(n) == (2 if n >= 1 else 0)
        assert no_edge.sequence.value(n) == (2 if n >= 2 else 0)


def test_kernel_and_noedge_exhaustion_budget_is_unknown():
    counter = ChainCounter(path_language(), joint_budget=1)
    p = make_graph({i: 'C' for i in range(4)}, [(0, 1, '-'), (2, 3, '-')])
    with pytest.raises(CountingUnknown, match='joint chain counting budget'):
        counter.count(p)
    # No partially computed count was cached as exact; independent requests
    # reset the operation budget and a connected query still works.
    assert counter.count(make_graph({0: 'C'})).sequence.value(10) == 11


def test_random_disconnected_patterns_against_injective_host_maps():
    rng = random.Random(170923)
    counter = ChainCounter(path_language())
    for trial in range(65):
        k = rng.randrange(2, 6)
        edges = [(i, i+1, '-') for i in range(k-1) if rng.random() < .35]
        p = make_graph({i: 'C' for i in range(k)}, edges)
        neg = tuple((u, v) for u in range(k) for v in range(u+1, k) if rng.random() < .15)
        adj = tuple(Adjacency(i, ('C',), ('-',), rng.choice(['=', '<', '>=']), rng.randrange(1, 4))
                    for i in range(k) if rng.random() < .25)
        result = counter.count(p, neg, adj)
        for n in [0, 1, 2, 4, 6, 8]:
            host = counter.language.instantiate(n)[0]
            assert result.sequence.value(n) == len(matches(p, host, neg, adj)), (trial, n, p, neg, adj)
