"""Finite host checks below are verification, not the synthesis reducer."""
from mfobs.core.graph import make_graph
from fractions import Fraction
from collections import Counter
from itertools import product
from pathlib import Path

import pytest

from mfobs.checks import eight_slots as reference
from mfobs.models import eight_slots as slot_counts
from mfobs.core.worklist import run
from mfobs.core.engine import mf
from mfobs.core.gml import rule_gml

from mfobs.inputs import InputError, load
from mfobs.core.linear_closure import RationalSpan
from mfobs.invariants.eight_slots import CONTRACT, CounterFromAxioms
from test_local_slot_preservation import count, identity, slot_host

ROOT = Path(__file__).resolve().parents[1]


def observable(word, cycle=False, extra_no_edges=()):
    n = len(word)
    edges = tuple((i, i+1) for i in range(n-1))
    if cycle:
        edges += ((n-1, 0),)
    graph = slot_host(tuple(int(x == 1) for x in word), edges)
    no_edges = tuple((3*i+1, 3*i+2) for i, x in enumerate(word) if x == 0)+extra_no_edges
    literal = identity(graph, no_edges)
    return literal, mf.Rule.fromGMLString(rule_gml(literal), add=False)


def sampled_hosts():
    for lengths in reference.TOPOLOGIES.values():
        edges = reference.topology_edges(lengths)
        for mask in (0, 1, 3, 15, 85, 129, 254, 255):
            bits = tuple((mask >> i) & 1 for i in range(8))
            yield lengths, bits, slot_host(bits, edges)


def test_possible_cycle_lengths_follow_from_size_and_degree_only():
    assert CounterFromAxioms(CONTRACT).lengths == (3, 4, 5, 8)
    with pytest.raises(ValueError):
        CounterFromAxioms({**CONTRACT, 'total_C': 9})
    with pytest.raises(ValueError):
        CounterFromAxioms({**CONTRACT, 'require_connected': True})


def test_single_slot_completion_requires_factor_one_half():
    _, rule = observable((1,))
    result = CounterFromAxioms(CONTRACT).reduce(rule)
    assert result.status == 'exact'
    assert len(result.polynomial) == 4 and set(result.polynomial.values()) == {Fraction(1, 2)}
    for lengths, bits, _ in sampled_hosts():
        assert slot_counts.polynomial_count(result.polynomial, bits, lengths) == sum(bits)


def test_free_is_not_unspecified_and_empty_count_uses_exact_eight():
    counter = CounterFromAxioms(CONTRACT)
    empty = counter.reduce(observable(())[1])
    free = counter.reduce(observable((0,))[1])
    any_site = counter.reduce(observable((None,))[1])
    assert free.polynomial != any_site.polynomial
    for lengths, bits, _ in sampled_hosts():
        assert slot_counts.polynomial_count(empty.polynomial, bits, lengths) == 1
        assert slot_counts.polynomial_count(free.polynomial, bits, lengths) == 8-sum(bits)
        assert slot_counts.polynomial_count(any_site.polynomial, bits, lengths) == 8


def test_completion_polynomials_preserve_raw_injective_counts():
    counter = CounterFromAxioms(CONTRACT)
    words = [word for n in range(1, 4) for word in product((None, 0, 1), repeat=n)]
    words += [tuple(state for _ in range(n)) for n in range(4, 9) for state in (None, 0, 1)]
    words += [tuple(i % 2 for i in range(n)) for n in range(4, 9)]
    checked = 0
    for word in words:
        for cycle in (False, True) if len(word) >= 3 else (False,):
            literal, rule = observable(word, cycle)
            result = counter.reduce(rule)
            assert result.status in ('exact', 'zero')
            for lengths, bits, host in sampled_hosts():
                direct = count(literal, host)
                assert direct == slot_counts.polynomial_count(result.polynomial, bits, lengths)
                checked += 1
    assert checked > 2000


@pytest.mark.parametrize('length', [6, 7, 9])
def test_impossible_cycle_component_has_a_structural_zero_certificate(length):
    result = CounterFromAxioms(CONTRACT).reduce(observable((None,)*length, cycle=True)[1])
    assert result.status == 'zero' and not result.polynomial


def test_noedge_between_path_ends_excludes_the_equal_length_completion():
    counter = CounterFromAxioms(CONTRACT)
    literal, rule = observable((0, None, 1), extra_no_edges=((0, 6),))
    result = counter.reduce(rule)
    assert {x['cycle_length'] for x in result.completions} == {4, 5, 8}
    for lengths, bits, host in sampled_hosts():
        assert count(literal, host) == slot_counts.polynomial_count(result.polynomial, bits, lengths)
    assert counter.reduce(observable((None,)*8, extra_no_edges=((0, 21),))[1]).status == 'zero'


def test_ninth_pattern_site_is_a_semantic_zero_not_an_exploration_budget():
    assert CounterFromAxioms(CONTRACT).reduce(observable((None,)*9)[1]).status == 'zero'


def test_disconnected_and_incomplete_patterns_remain_unknown():
    counter = CounterFromAxioms(CONTRACT)
    for graph in (slot_host((0, 0)), make_graph({0: 'C'})):
        result = counter.reduce(mf.Rule.fromGMLString(rule_gml(identity(graph)), add=False))
        assert result.status == 'unknown'


def test_sparse_elimination_tracks_actual_pattern_coordinates():
    span = RationalSpan()
    a = {'a': Fraction(2), 'b': Fraction(3)}
    b = {'b': Fraction(2), 'c': Fraction(-1)}
    assert span.insert('P0', a)
    assert span.insert('P1', b)
    rest, coefficients = span.reduce({'a': Fraction(6), 'b': Fraction(5), 'c': Fraction(2)})
    assert not rest and coefficients == {'P0': 3, 'P1': -2}
    assert not span.insert('duplicate', a)


def test_actual_discovery_is_independent_of_host_or_orbit_evaluation(monkeypatch):
    from mfobs.models import ring_language
    from mfobs.core import mixtures
    from mfobs.checks import local_slots_reference
    def forbidden(*args, **kwargs):
        raise AssertionError('Synthesis called a host/state oracle')
    for name in ('orbit_generator', 'slot_ring_graph', 'evaluate_observable', 'word_values'):
        monkeypatch.setattr(ring_language, name, forbidden)
    monkeypatch.setattr(mixtures, 'fast_count', forbidden)
    monkeypatch.setattr(local_slots_reference, 'typed_match_count', forbidden)
    monkeypatch.setattr(slot_counts, 'polynomial_count', forbidden)
    inputs = load(ROOT/'examples/eight-slots/model.json')
    result = run(inputs, max_expansions=3)
    assert result['basis'] == ['P0', 'P1', 'P2', 'P3']
    assert result['equations_by_reaction']['P0'] == {
        'basal_bind': {'P1': '1'}, 'cooperative_bind': {'P2': '1'}, 'unbind': {'P0': '-1'}}
    monkeypatch.undo()  # Independent verification starts after synthesis.
    for lengths, bits, host in sampled_hosts():
        counts = {n: slot_counts.polynomial_count(v, bits, lengths)
                  for n, v in result['basis_polynomials'].items()}
        occupied = sum(bits)
        neighbours = sum(bits[a] != bits[b] for a, b in reference.topology_edges(lengths))
        expected = {'k_on': 8 - occupied, 'k_coop': neighbours, 'k_off': -occupied}
        for channel, row in result['equations_by_channel']['P0'].items():
            assert sum(Fraction(c) * counts[n] for n, c in row.items()) == expected[channel]
    assert not result['closed'] and result['matrix'] is None
    assert result['synthesis_state_evaluations'] == result['synthesis_orbit_evaluations'] == 0
    assert result['initial_values'] == {}


@pytest.mark.parametrize('policy', [
    {'connected': True}, {'include_empty': False}, {'mode': 'restricted'},
])
def test_uncertified_overlap_shortcuts_are_rejected(policy):
    inputs = load(ROOT/'examples/eight-slots/model.json')
    inputs.data['composition'].update(policy)
    with pytest.raises(InputError):
        run(inputs, max_expansions=1)


@pytest.mark.parametrize('reaction', ['basal_bind', 'cooperative_bind', 'unbind'])
def test_omega_bounds_prune_compositions_without_changing_generator(reaction):
    from mfobs.core.engine import commutator
    inputs = load(ROOT/'examples/eight-slots/model.json')
    adapter = inputs.adapter()
    query = observable((1,))[1]
    rule = inputs.engine_rules()[reaction]
    full = commutator(query, rule)
    pruned = commutator(query, rule, constraints=adapter.composition_constraints)

    def polynomial(image):
        result = Counter()
        for key, coefficient in image.observable_coefficients.items():
            reduced = adapter.counter.reduce(image.observable_representatives[key])
            assert reduced.status in ('exact', 'zero')
            for monomial, value in reduced.polynomial.items():
                result[monomial] += coefficient * value
        return {key: value for key, value in result.items() if value}

    assert polynomial(pruned) == polynomial(full)
    assert pruned.left_raw + pruned.right_raw < full.left_raw + full.right_raw
