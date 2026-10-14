"""Structural, generator, and report checks for the NoEdge oscillator."""
from mfobs.core.molecules import BOND_ORDERS as BOND_ORDER
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import edge_rows, make_graph, node_labels, node_rows, same_graph
from mfobs.core.rules import added_edges, rule_graph
from mfobs.core.conditions import engine_no_edges

import mfobs.checks.growing_ring as ring_checks
from itertools import product
from pathlib import Path

import pytest

from mfobs.core.conditions import matches
from mfobs.core.graph import apply_match
from mfobs.models import structural_nac_oscillator as chemistry
from mfobs.models import structural_nac_words as words


ROOT = Path(__file__).resolve().parents[1]


def bond_sums(graph):
    result = {node: 0 for node in node_labels(graph)}
    for left, right, label in edge_rows(graph):
        result[left] += BOND_ORDER[label]
        result[right] += BOND_ORDER[label]
    return result


def test_unbounded_family_uses_only_ordinary_chemical_labels_and_full_valence():
    for states in ((0, 1, 2), (0,) * 12, (0, 2, 1, 1, 0, 2, 2, 1)):
        graph = chemistry.macrocycle(states)
        assert set(node_labels(graph).values()) == {"C", "H", "O", "S"}
        assert {label for _, _, label in edge_rows(graph)} <= {"-", "="}
        assert all(
            bond_sums(graph)[node] == chemistry.VALENCE[label]
            for node, label in node_rows(graph)
        )
    assert chemistry.macrocycle((0,) * 13).numVertices > chemistry.macrocycle((0,) * 12).numVertices
    for state in range(3):
        fragment = chemistry.capped_site(state)
        assert all(
            bond_sums(fragment)[node] == chemistry.VALENCE[label]
            for node, label in node_rows(fragment)
        )
        assert to_rdkit(fragment).GetNumAtoms() == fragment.numVertices


def test_all_dpo_rules_are_atom_and_local_bond_order_balanced():
    from mfobs.core.engine import mf
    from mfobs.core.gml import rule_gml

    assert len(chemistry.rules()) == 15
    for name, rule in chemistry.rules().items():
        assert chemistry.atom_balance(rule), name
        assert chemistry.local_bond_order_balance(rule), name
        assert node_rows(rule_graph(rule, "left")) == node_rows(rule_graph(rule, "right"))
        assert {label for _, _, label in (*edge_rows(rule_graph(rule, "left")), *edge_rows(rule_graph(rule, "right")))} <= {"-", "="}
        identity = {node: node for node in node_labels(rule_graph(rule, "left"))}
        assert same_graph(apply_match(rule, rule_graph(rule, "left"), identity), rule_graph(rule, "right"))
        assert mf.Rule.fromGMLString(rule_gml(rule), add=False) is not None


def test_per_match_factors_exactly_remove_mapped_hydrogen_choices():
    for name, rule in chemistry.rules().items():
        complete_host = chemistry.complete_with_h(rule_graph(rule, "left"))
        raw = len(matches(rule, complete_host))
        assert chemistry.rate_factors()[name] * raw == 1


def test_each_noedge_is_active_on_an_unchanged_probe():
    guarded = [rule for name, rule in chemistry.rules().items()
               if name in chemistry.RULE_TRANSITIONS and chemistry.RULE_TRANSITIONS[name][3] is not None]
    assert len(guarded) == 3
    for rule in guarded:
        assert len(matches(rule, rule_graph(rule, "left"))) == 1
        product_edges = {(u, v) for u, v, _ in added_edges(rule)}
        guards = [pair for pair in engine_no_edges(rule) if pair not in product_edges]
        assert len(guards) == 1
        left, right = guards[0]
        blocked = make_graph(node_labels(rule_graph(rule, "left")), (*edge_rows(rule_graph(rule, "left")), (left, right, "-")))
        assert not matches(rule, blocked)
        # The guarded probe is context, not the bond being changed.
        assert (left, right, "-") not in edge_rows(rule_graph(rule, "right"))


def test_swap_family_preserves_every_one_site_observable():
    for reaction in chemistry.SWAP_REACTIONS:
        assert all(words.generator_row((state,), reaction) == {} for state in range(3))


def test_growth_copies_template_changes_size_and_cancels_from_one_site_fractions():
    state = (0, 1, 2, 0, 2, 1, 0)
    for reaction, copied in chemistry.GROW_STATES.items():
        events = list(chemistry.transition_events(state, reaction))
        assert len(events) == state.count(copied)
        assert all(len(target) == len(state) + 1 for target in events)
        assert all(target.count(copied) == state.count(copied) + 1 for target in events)
        for observed in range(3):
            row = words.generator_row((observed,), reaction)
            assert row == ({(copied,): 1} if observed == copied else {})
    # Summing the three beta channels gives n_i'=beta*n_i and N'=beta*N,
    # hence zero contribution to d(n_i/N)/dt.
    for observed in range(3):
        raw = sum(
            ring_checks.evaluate_row(words.generator_row((observed,), reaction), state)
            for reaction in chemistry.GROW_REACTIONS
        )
        assert raw == chemistry.word_count(state, (observed,))


def test_sparse_growth_row_matches_exhaustive_context_enumeration():
    from collections import Counter

    def exhaustive(word, reaction):
        copied = chemistry.GROW_STATES[reaction]
        if len(word) == 1:
            return {(copied,): 1} if word == (copied,) else {}
        length = len(word)
        context_length = 2 * length - 2
        anchor = length - 2
        free = [index for index in range(context_length) if index != anchor]
        result = Counter()
        for values in product(range(3), repeat=len(free)):
            context = [None] * context_length
            context[anchor] = copied
            for index, value in zip(free, values):
                context[index] = value
            before = tuple(context)
            after = before[:anchor + 1] + (copied,) + before[anchor + 1:]
            coefficient = sum(
                after[start:start + length] == word for start in range(length)
            ) - sum(
                before[start:start + length] == word
                for start in range(length - 1)
            )
            if coefficient:
                result[before] += coefficient
        return {key: value for key, value in sorted(result.items()) if value}

    for word in ((0, 1), (1, 2, 0), (2, 0, 1, 2), (0, 2, 1, 0, 2)):
        for reaction in chemistry.GROW_REACTIONS:
            assert words.growth_generator_row(word, reaction) == exhaustive(
                word, reaction
            )


def test_literal_growth_rule_consumes_complete_precursor_and_releases_water():
    def shifted(graph, offset):
        return make_graph(
            {node + offset: label for node, label in node_rows(graph)},
            [(left + offset, right + offset, bond) for left, right, bond in edge_rows(graph)],
        )

    state = (0, 1, 2, 0, 2, 1, 0)
    ring = chemistry.macrocycle(state)
    precursor = shifted(chemistry.capped_site(0), 1000)
    host = make_graph({**node_labels(ring), **node_labels(precursor)}, (*edge_rows(ring), *edge_rows(precursor)))
    rule = chemistry.rules()[chemistry.GROW_0]
    growth_matches = matches(rule, host)
    assert len(growth_matches) == state.count(0)
    product = apply_match(rule, host, growth_matches[0])

    adjacency = {node: set() for node in node_labels(product)}
    for left, right, _bond in edge_rows(product):
        adjacency[left].add(right)
        adjacency[right].add(left)
    components = []
    unseen = set(adjacency)
    while unseen:
        component, frontier = set(), [next(iter(unseen))]
        while frontier:
            node = frontier.pop()
            if node in component:
                continue
            component.add(node)
            frontier.extend(adjacency[node] - component)
        unseen -= component
        components.append(component)
    water = next(component for component in components if len(component) == 3)
    assert sorted(node_labels(product)[node] for node in water) == ["H", "H", "O"]
    enlarged = next(component for component in components if len(component) != 3)
    assert len(enlarged) == chemistry.macrocycle(state + (0,)).numVertices


@pytest.mark.parametrize("ring_size", [5, 6])
def test_exact_word_rows_equal_direct_event_sums(ring_size):
    for states in product(range(3), repeat=ring_size):
        for observable in (*product(range(3), repeat=1), *product(range(3), repeat=2)):
            for reaction in chemistry.REACTION_ORDER:
                row = words.generator_row(observable, reaction)
                assert ring_checks.evaluate_row(row, states) == ring_checks.direct_jump_sum(
                    observable, reaction, states
                )


def test_noedge_generates_triplet_observables_and_mixing_generates_pairs():
    assert words.generator_row((0,), chemistry.GUARD_10) == {
        (1, 0, 0): -1,
        (1, 0, 2): -1,
    }
    assert words.generator_row((0, 1), chemistry.SWAP_01)
    assert chemistry.one_site_rhs((1 / 3, 1 / 3, 1 / 3), 1, .12) == pytest.approx(
        (0, 0, 0)
    )
    eigenvalues = chemistry.one_site_eigenvalues(1, .12)
    assert eigenvalues[0].real == pytest.approx(-.02)
    assert eigenvalues[0].imag == pytest.approx(.611991, rel=1e-5)


def test_budgeted_worklist_discovers_connected_graph_words_and_is_unknown():
    from mfobs.core.worklist import run
    from mfobs.inputs import load

    inputs = load(ROOT / "examples" / "structural-nac-oscillator" / "model.json")
    result = run(inputs, max_expansions=6)
    assert result["status"] == "UNKNOWN"
    assert result["basis"][:3] == ["Q0", "Q1", "Q2"]
    from fractions import Fraction
    forms = {n: {tuple(map(int, w[1:])): Fraction(c) for w, c in v['terms'].items()}
             for n, v in result['word_coordinates'].items()}
    assert any(len(v) > 1 for v in forms.values())
    # Check the new *combined* coordinates on actual cyclic-word transitions.
    # generator_row already counts normalized word events; mapped hydrogen
    # multiplicities have been removed in that theorem, not a second time here.
    for host in [(0, 1, 2) * 4, (0, 0, 1, 2, 1, 1, 0, 2, 2, 0, 1, 2)]:
        value = lambda form, word: sum(c * chemistry.word_count(word, w) for w, c in form.items())
        for step in result['algorithm1_steps']:
            if 'coordinates' not in step:
                continue
            reactions = [r for r in inputs.data['rules'] if r['rate'] == step['channel']]
            direct = sum(value(forms[step['source']], target) - value(forms[step['source']], host)
                         for reaction in reactions
                         for target in chemistry.transition_events(host, reaction['id']))
            expected = sum(Fraction(c) * value(forms[n], host) for n, c in step['coordinates'].items())
            assert direct == expected, (step['source'], step['channel'])
    assert result['obligations'][0]['source'] in result['pending']
    assert result["raw_compositions"] == 0
    assert result["semantic_contract"]["all_rules_verified"]


