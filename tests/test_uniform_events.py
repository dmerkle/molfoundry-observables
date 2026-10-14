"""Event extraction must be independent of terminal paths and fail closed."""
from types import SimpleNamespace

import sympy as sp

from mfobs.core.pattern_support import Pattern, Support
from mfobs.core.uniform_events import EventClass, extract_uniform_events, constant_partition


def extract(image, activity=1, partition=constant_partition, **kwargs):
    support = Support('POLYNOMIAL_SUPPORT_CLOSED',
        {'x': Pattern('x', None, 'arbitrary-count')},
        {'x': {'r': sp.sympify(image)}}, {'r': sp.sympify(activity)}, {}, [], [])
    return extract_uniform_events(support, [SimpleNamespace(identifier='r', rate='k')],
        partition=partition, scope='test family',
        uniform_match_argument=kwargs.get('argument', 'declared uniform matches'),
        partition_argument='declared exhaustive enabled-state partition')


def test_nonboolean_partition_and_rate_reconstruction():
    # A backend with three possible count values; no suffix/gate conventions.
    x = sp.Symbol('x')
    def ternary(source, jumps):
        return [EventClass({x: i}, sp.prod((x-j)/(i-j) for j in range(3) if j != i),
                           {'x': i}) for i in range(3)]
    result = extract(3*(x+1), 3, ternary)
    assert result['status'] == 'FINITE_EVENTS'
    assert [e['jump']['x'] for e in result['events']] == [1, 2, 3]
    rates = [sp.sympify(e['activity']['text']) for e in result['events']]
    assert sp.expand(sum(rates)) == 3
    assert sp.expand(sum(e['jump']['x']*a for e,a in zip(result['events'],rates))) == 3*(x+1)
    for i in range(3):
        assert [a.subs(x,i) for a in rates] == [3 if j==i else 0 for j in range(3)]


def test_missing_justification_cannot_certify_events():
    assert extract(1, argument='')['status'] == 'UNKNOWN_UNIFORM_MATCH_OBLIGATION'


def test_nonconstant_or_nonintegral_jump_cannot_certify_events():
    for image in (sp.Symbol('x'), sp.Rational(1,2)):
        result = extract(image)
        assert result['status'] == 'UNKNOWN_EVENT_EXTRACTION'
        assert result['events'] == []


def test_empty_partition_is_not_a_closed_empty_event_model():
    assert extract(1, partition=lambda source,jumps: [])['status'] == 'UNKNOWN_EVENT_EXTRACTION'


def test_zero_activity_does_not_divide_by_zero():
    result = extract(0, 0)
    assert result['status'] == 'FINITE_EVENTS'
    assert result['events'] == []
