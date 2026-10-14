"""Scientific regressions distinguishing Algorithm 1 from pattern saturation."""
from fractions import Fraction as F
from unittest.mock import patch

import pytest

from molfoundry import graphGMLString, ruleGMLString

from mfobs import Observable, RatedRule, SimpleGraphs, SimpleGraphBackend, ModelInputs, analyze
from mfobs.core.jump_closure_worklist import RuleExpression, discover_observables_by_channel
from mfobs.core.linear_closure import close_linear_span
from mfobs.core.polynomial_jump_closure import (
    JumpEvent, JumpModelProvenance, PolynomialJumpModel, PolynomialRing,
    discover_linear_observables_from_events as discover_linear_observables, generator_action,
)


def model(ring, events):
    return PolynomialJumpModel('regression', ring, tuple(events),
        {'status': 'VERIFIED_ALL_STATES'}, 'regression',
        JumpModelProvenance('test', 'explicit CTMC', 'test deltas', 'test rates',
                            'none', 'none', False, True))


def event(name, channel, delta, rate):
    return JumpEvent(name, channel, delta, rate, name, 'enabled by propensity')


def test_whole_residual_closes_where_monomial_saturation_does_not():
    ring = PolynomialRing(('q', 'x', 'y'))
    q, x, y = map(ring.coordinate, ring.variables)
    process = model(ring, [event('q_x', 'k', (1, 0, 0), x),
                           event('q_y', 'k', (1, 0, 0), y),
                           event('convert', 'k', (0, 1, -1), ring.multiply(x, y))])
    # A valid CTMC on nonnegative counts: conversion is disabled at y=0.
    # Splitting x+y exposes the unbounded x^a*y moment hierarchy.
    result = discover_linear_observables(process, ['q'], maximum_dimension=2)
    assert result['status'] == 'CLOSED'
    assert result['dimension'] == 2
    assert generator_action(process, q) == ring.add(x, y)
    assert generator_action(process, ring.add(x, y)).is_zero
    assert result['steps'][0]['coordinates'] == {result['basis'][1]: '1'}
    assert result['steps'][1]['coordinates'] == {}


def test_same_rate_rules_cancel_before_span_extension():
    ring = PolynomialRing(('x', 'y'))
    propensity = ring.multiply(*(ring.coordinate(n) for n in ring.variables))
    process = model(ring, [event('birth', 'k', (1, 0), propensity),
                           event('death', 'k', (-1, 0), propensity)])
    assert discover_linear_observables(process, ['x'], maximum_dimension=1)['status'] == 'CLOSED'
    separate = model(ring, [event('birth', 'a', (1, 0), propensity),
                            event('death', 'b', (-1, 0), propensity)])
    assert discover_linear_observables(separate, ['x'], maximum_dimension=1)['status'] == 'UNKNOWN_DIMENSION_LIMIT'


def test_budget_retains_failed_channel_image_and_entire_frontier():
    result = close_linear_span(seeds={'x': {'x': F(1)}, 'y': {'y': F(1)}},
        channels=['a', 'b'], image=lambda n, c: ({'z': F(1)}, {}), maximum_dimension=2)
    assert result.pending == ['x', 'y']
    assert [(o['source'], o['channel']) for o in result.obligations] == [
        ('x', 'a'), ('x', 'b'), ('y', 'a'), ('y', 'b')]
    assert result.obligations[0]['residual'] == {'z': '1'}
    assert result.steps[0]['status'] == 'UNKNOWN_DIMENSION_LIMIT'


def birth_input():
    omega = SimpleGraphBackend(SimpleGraphs(vertex_labels={'A'}, edge_labels={'-'}))
    birth = RatedRule(ruleGMLString(
        'rule [ ruleID "birth" right [ node [ id 0 label "A" ] ] ]', add=False), 'k')
    query = Observable('A', graphGMLString(
        'graph [ node [ id 0 label "A" ] ]', add=False))
    return omega, birth, query


def test_real_graph_composition_preserves_scaled_query():
    omega, birth, query = birth_input()
    seed = omega._engine_observable(query)
    result = discover_observables_by_channel(
        seeds={'twice_A': RuleExpression.from_rules([seed]).scale(F(2))},
        generators={'k': RuleExpression.from_rules([birth.rule])},
        pair_commutator=omega._pair, maximum_dimension=2,
        jump_closure=__import__('mfobs.core.engine', fromlist=['dpo_jump_closure']).dpo_jump_closure)
    assert result.status == 'CLOSED'
    assert set(result.basis['twice_A'].coefficients.values()) == {F(2)}
    assert result.steps[0]['coordinates'] == {'O1': '2'}


def test_public_api_keeps_independent_branch_results(capsys):
    omega, birth, query = birth_input()
    result = analyze(ModelInputs((omega).omega, {r.identifier:r for r in iter([birth])}, {(query).name:query}), backend=omega, maximum_dimension=1, uniform_matches=True)
    assert result.linear.status == 'UNKNOWN_DIMENSION_LIMIT'
    assert result.linear.data["pending"] == ['A']
    assert result.support.status == 'POLYNOMIAL_SUPPORT_CLOSED'
    assert result.pjc.status == 'CLOSED'
    result.show()
    assert 'UNKNOWN_DIMENSION_LIMIT' in capsys.readouterr().out
    missing = analyze(ModelInputs((omega).omega, {r.identifier:r for r in [birth]}, {(query).name:query}), backend=omega)
    assert missing.linear.status == 'CLOSED'
    assert missing.pjc.status == 'UNKNOWN_UNIFORMITY'
    assert not missing.record()["events"]
    missing.show()


def test_public_algorithm_calls_run_only_the_selected_search():
    from mfobs.backends import simple_graphs

    omega, birth, query = birth_input()
    with (
        patch.object(simple_graphs, 'discover_observables_by_channel',
                     wraps=simple_graphs.discover_observables_by_channel) as linear_search,
        patch.object(simple_graphs, 'discover_support',
                     wraps=simple_graphs.discover_support) as support_search,
        patch('mfobs.core.polynomial_jump_closure.discover_linear_observables_from_events',
              wraps=discover_linear_observables) as event_linear_search,
    ):
        jump = omega.jump_closure(ModelInputs((omega).omega, {r.identifier:r for r in [birth]}, {(query).name:query}), uniform_matches=True)
        assert jump.pjc['status'] == 'CLOSED'
        assert support_search.call_count == 1
        assert linear_search.call_count == event_linear_search.call_count == 0

        linear = omega.linear_closure(ModelInputs((omega).omega, {r.identifier:r for r in [birth]}, {(query).name:query}), maximum_dimension=1)
        assert linear.status == 'UNKNOWN_DIMENSION_LIMIT'
        assert linear_search.call_count == 1
        assert support_search.call_count == 1
        assert event_linear_search.call_count == 0


@pytest.mark.parametrize('limits, status', [
    ({'maximum_support_patterns': 1}, 'UNKNOWN_PATTERN_LIMIT'),
    ({'maximum_support_expansions': 0}, 'UNKNOWN_EXPANSION_LIMIT'),
])
def test_standalone_jump_closure_retains_unfinished_support(limits, status):
    omega, birth, query = birth_input()
    result = omega.jump_closure(ModelInputs((omega).omega, {r.identifier:r for r in [birth]}, {(query).name:query}), uniform_matches=True, **limits)
    assert result.support.status == result.pjc["status"] == status
    assert result.support.pending
    assert result.compilation is None


def test_counting_adapter_sums_normalized_rule_images_before_extension():
    from mfobs.core.engine import mf
    from mfobs.core.identity import Adapter
    from mfobs.core.worklist import State
    from types import SimpleNamespace as NS

    seed = mf.Rule.fromGMLString('rule [ context [ node [ id 0 label "X" ] ] ]', add=False)
    xy = mf.Rule.fromGMLString('rule [ context [ node [ id 0 label "X" ] node [ id 1 label "Y" ] ] ]', add=False)
    reactions = [{'id': 'birth', 'rate': 'k', 'per_match': '1/2'},
                 {'id': 'death', 'rate': 'k', 'per_match': '1'}]
    inputs = NS(data={'initial_observables': ['X'], 'rules': reactions},
                engine_observables=lambda: {'X': seed})
    adapter = Adapter(inputs)
    state = State(adapter)
    class Composer:
        records = []
        def pair(self, source, rule, reaction):
            # Exact images of birth/death at XY; first raw match count is 2XY.
            record = {'source': source, 'reaction': reaction['id']}
            self.records.append(record)
            coefficient = 2 if reaction['id'] == 'birth' else -1
            return NS(observable_coefficients={xy.canonKey: coefficient},
                      observable_representatives={xy.canonKey: xy}), record
    composer = Composer()
    adapter.initialize(state, composer)
    assert adapter.run_closure(state, composer, max_expansions=1, maximum_dimension=1)
    assert state.channel_equations == {'X': {'k': {}}}
    assert state.basis == ['X']
    # Per-rule images are nonzero and outside the span; they must not be
    # exported as zero rows merely because their channel sum vanishes.
    assert state.equations == {}


def test_simple_representative_changes_coordinates_but_not_the_extended_span():
    result = close_linear_span(seeds={'x': {'x': F(1)}}, channels=['k'],
        image=lambda n, c: ({'x': F(3), 'y': F(1)}, {}),
        normalize=lambda residual: ({'x': F(1), 'y': F(1)}, F(1)),
        maximum_expansions=1)
    assert result.steps[0]['coordinates'] == {'x': '2', 'O1': '1'}
    assert result.steps[0]['representative_coordinates'] == {'x': '1'}
    assert result.basis['O1'] == {'x': F(1), 'y': F(1)}
    assert result.expansions[0]['before'] == ['x']
    assert result.expansions[0]['after'] == ['x', 'O1']


def test_extension_cannot_smuggle_in_an_unrelated_direction():
    with pytest.raises(AssertionError, match='exactly the residual direction'):
        close_linear_span(seeds={'x': {'x': F(1)}}, channels=['k'],
            image=lambda n, c: ({'y': F(1), 'z': F(1)}, {}),
            normalize=lambda residual: ({'y': F(1)}, F(1)))
