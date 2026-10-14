"""Independent adversarial checks of the joint-counter extension.

Use a labelled, edge-labelled, nonuniformly decorated grammar with a special
zero graph; the host oracle enumerates complete injective maps independently.
"""
from mfobs.core.graph import make_graph, edge_rows, has_edge, node_labels
from itertools import product
import random



from mfobs.core.conditions import Adjacency, matches
from mfobs.invariants.chain import ChainLanguage, ChainCounter, components


def decorated_language():
    # The same C port receives both joins in B. Boundary decorations differ,
    # and zero has an edge label not appearing at any positive length.
    return ChainLanguage({'kind': 'one-port-decorated-path', 'join_bond': '-',
        'zero': {'nodes': [[0, 'C'], [1, 'N'], [2, 'C']],
                 'edges': [[0, 1, 'z'], [1, 2, '=']]},
        'left': {'nodes': [[0, 'C'], [1, 'O']], 'edges': [[0, 1, '=']],
                 'incoming': None, 'outgoing': 0},
        'repeat': {'nodes': [[0, 'C'], [1, 'O'], [2, 'N']],
                   'edges': [[0, 1, '-'], [1, 2, '=']],
                   'incoming': 0, 'outgoing': 0},
        'right': {'nodes': [[0, 'C'], [1, 'O'], [2, 'N']],
                  'edges': [[0, 1, '-'], [0, 2, 'e']],
                  'incoming': 0, 'outgoing': None}})


def assert_oracle(counter, pattern, no_edges=(), adjacency=(), lengths=(0, 1, 2, 3, 5, 8)):
    result = counter.count(pattern, no_edges, adjacency)
    for n in lengths:
        host, _ = counter.language.instantiate(n)
        want = len(matches(pattern, host, no_edges, adjacency))
        assert result.sequence.value(n) == want, (n, pattern, no_edges, adjacency, result.witness)
    return result


def test_multiple_edge_signatures_and_conflicting_quotient_edges():
    counter = ChainCounter(decorated_language())
    # Same endpoints' labels but different positive edge labels. A quotient
    # identifying both ends would require parallel distinct edges and is zero.
    pattern = make_graph({0: 'C', 1: 'O', 2: 'C', 3: 'O'},
                         [(0, 1, '='), (2, 3, '-')])
    result = assert_oracle(counter, pattern)
    assert result.sequence.value(0) == 0
    assert_oracle(counter, pattern, ((0, 2), (1, 3)))


def test_true_kernel_enumeration_with_conflicting_parallel_edge_quotients():
    # Symmetric alternating square gives each C and O both edge labels. Thus
    # supports overlap on BOTH endpoints. Kernel candidates identifying both
    # ends must actually reach the incompatible-parallel-edge rejection.
    square = {'nodes': [[0, 'C'], [1, 'O'], [2, 'C'], [3, 'O']],
              'edges': [[0, 1, '-'], [1, 2, '='], [2, 3, '-'], [3, 0, '=']]}
    language = ChainLanguage({'kind': 'one-port-decorated-path', 'join_bond': '-', 'zero': square,
        'left': {**square, 'incoming': None, 'outgoing': 0},
        'repeat': {**square, 'incoming': 0, 'outgoing': 0},
        'right': {**square, 'incoming': 0, 'outgoing': None}})
    counter = ChainCounter(language)
    pattern = make_graph({0: 'C', 1: 'O', 2: 'C', 3: 'O'},
                         [(0, 1, '='), (2, 3, '-')])
    result = assert_oracle(counter, pattern)
    assert result.witness['method'] == 'component-kernel-partition'
    # In one square, every = edge meets every - edge, hence zero injective
    # disjoint pair; positive lengths admit pairs from different tiles.
    assert result.sequence.value(0) == 0
    for n in (1, 2, 1000):
        assert result.sequence.value(n) == 4*n*(n+1)


def test_zero_only_edge_label_in_cross_noedge_complement():
    counter = ChainCounter(decorated_language())
    pattern = make_graph({0: 'C', 1: 'N'})
    result = assert_oracle(counter, pattern, ((0, 1),))
    assert result.witness['method'] == 'cross-noedge-complement'
    assert 'z' in {part['bond'] for part in result.witness['edge_corrections']}
    assert result.sequence.value(0) == 0


def test_colliding_centres_keep_conjunction_of_different_adjacency_predicates():
    counter = ChainCounter(decorated_language())
    pattern = make_graph({0: 'C', 1: 'C', 2: 'C'})
    predicates = (Adjacency(0, ('O',), ('-', '='), '>=', 1),
                  Adjacency(1, ('C',), ('-',), '<=', 1),
                  Adjacency(2, ('N',), ('e', 'z'), '=', 0))
    result = assert_oracle(counter, pattern, (), predicates)
    assert result.witness['method'] == 'component-kernel-partition'
    assert_oracle(counter, pattern, ((0, 1), (1, 2)), predicates)


def test_inconsistent_colliding_predicates_contribute_zero_not_negative_count():
    counter = ChainCounter(decorated_language())
    pattern = make_graph({0: 'C', 1: 'C'})
    predicates = (Adjacency(0, ('C',), ('-',), '=', 0),
                  Adjacency(1, ('C',), ('-',), '>=', 1))
    assert_oracle(counter, pattern, (), predicates)


def test_free_vertices_mixed_labels_conditioned_core_and_special_zero():
    counter = ChainCounter(decorated_language())
    pattern = make_graph({0: 'C', 1: 'O', 2: 'C', 3: 'N', 4: 'N'}, [(0, 1, '-')])
    predicates = (Adjacency(0, ('C',), ('-',), '=', 2),)
    result = assert_oracle(counter, pattern, (), predicates)
    assert result.witness['method'] == 'free-vertex-falling-factorials'
    assert_oracle(counter, make_graph({0: 'C', 1: 'N', 2: 'C'}), ((0, 1),))


def test_all_two_vertex_label_noedge_and_unary_combinations():
    counter = ChainCounter(decorated_language())
    for labels in product(('C', 'O', 'N', 'ABSENT'), repeat=2):
        pattern = make_graph(dict(enumerate(labels)))
        for no_edge in [(), ((0, 1),)]:
            for predicates in [(), (Adjacency(0, ('C', 'N'), ('-', '=', 'e', 'z'), '<', 2),)]:
                assert_oracle(counter, pattern, no_edge, predicates)


def test_random_mixed_label_disconnected_patterns_with_external_neighbour_conditions():
    rng = random.Random(593108)
    counter = ChainCounter(decorated_language())
    checked = 0
    for trial in range(140):
        k = rng.randrange(2, 6)
        nodes = {2*i+7: rng.choice(('C', 'O', 'N')) for i in range(k)}
        ids = list(nodes)
        edges = [(u, v, rng.choice(('-', '=', 'e', 'z')))
                 for i, u in enumerate(ids) for v in ids[i+1:]
                 if rng.random() < .17]
        pattern = make_graph(nodes, edges)
        if len(components(pattern)) == 1:
            continue
        negative = tuple((u, v) for i, u in enumerate(ids) for v in ids[i+1:]
                         if rng.random() < .2)
        predicates = tuple(Adjacency(u, rng.choice([('C',), ('O',), ('N',), ('C', 'N')]),
                                    rng.choice([('-',), ('=',), ('e', 'z'), ('-', '=')]),
                                    rng.choice(['=', '<', '<=', '>', '>=']), rng.randrange(4))
                           for u in ids if rng.random() < .35)
        assert_oracle(counter, pattern, negative, predicates, lengths=(0, 1, 2, 4, 7))
        checked += 1
    assert checked >= 100


def test_random_positive_witness_patterns_with_admissible_cross_guards():
    rng = random.Random(761204)
    counter = ChainCounter(decorated_language())
    checked = 0
    for trial in range(80):
        source, _ = counter.language.instantiate(rng.randrange(2, 6))
        ids = sorted(rng.sample(list(node_labels(source)), rng.randrange(2, 6)))
        edges = [(u, v, b) for u, v, b in edge_rows(source)
                 if u in ids and v in ids and rng.random() < .7]
        pattern = make_graph({u: node_labels(source)[u] for u in ids}, edges)
        if len(components(pattern)) == 1:
            continue
        negative = tuple((u, v) for i, u in enumerate(ids) for v in ids[i+1:]
                         if not has_edge(source, u, v) and rng.random() < .3)
        predicates = []
        for u in ids:
            if rng.random() < .3:
                nodes = rng.choice([('C',), ('N', 'O')])
                bonds = rng.choice([('-',), ('=', 'e', 'z')])
                probe = Adjacency(u, nodes, bonds, '=', 0)
                predicates.append(Adjacency(u, nodes, bonds, '=', probe.count_at(source, u)))
        result = assert_oracle(counter, pattern, negative, tuple(predicates), lengths=(0, 1, 2, 4, 7))
        assert result.sequence, 'Patterns were extracted with an admissible source embedding'
        checked += 1
    assert checked >= 60
