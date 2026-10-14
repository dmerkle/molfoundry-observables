"""Invariant graph queries give a certified zero generator and constant dynamics."""
import numpy as np
import pytest
import molfoundry as mf

from mfobs import ModelInputs, Observable, RatedRule, SimpleGraphs, SimpleGraphBackend
from mfobs.core.polynomial_jump_closure import (
    discover_linear_observables_from_events, generator_action, one_step_event_expansion,
)


@pytest.fixture(scope="module")
def invariant_query():
    rule = mf.Rule.fromGMLString('''rule [ ruleID "subdivide"
        context [node[id 0 label "L"] node[id 1 label "R"]]
        left [edge[source 0 target 1 label "-"]]
        right [node[id 2 label "L"] edge[source 0 target 2 label "-"]
               edge[source 2 target 1 label "-"]] ]''', add=False)
    query = mf.Graph.fromGMLString('''graph [node[id 0 label "L"] node[id 1 label "R"]
        edge[source 0 target 1 label "-"]]''', add=False)
    omega = SimpleGraphs(vertex_labels={"L", "R"}, edge_labels={"-"})
    return SimpleGraphBackend(omega), ModelInputs(omega, {rule.name: RatedRule(rule)},
                                                 {"E_LR": Observable("E_LR", query)})


def test_invariant_graph_query_closes_both_algorithms(invariant_query):
    backend, inputs = invariant_query
    linear = backend.linear_closure(inputs)
    assert linear.status == "CLOSED"
    assert linear.steps[0]["coordinates"] == {}
    jump = backend.jump_closure(inputs, uniform_matches=True)
    assert jump.support.status == "POLYNOMIAL_SUPPORT_CLOSED"
    assert jump.compilation.model.events == ()
    assert jump.pjc["status"] == "CLOSED"
    assert jump.pjc["operator_span"]["dimension"] == 0
    assert jump.pjc["operator"]["terms"] == []
    assert jump.pjc["first_moment_generator_images"] == {"E_LR": {}}
    assert backend.jump_closure(inputs).pjc["status"] == "UNKNOWN_UNIFORMITY"


def test_zero_event_generator_and_linear_closure(invariant_query):
    backend, inputs = invariant_query
    model = backend.jump_closure(inputs, uniform_matches=True).compilation.model
    assert generator_action(model, model.ring.coordinate("E_LR")).is_zero
    assert one_step_event_expansion(model, {"E_LR": 3})["transitions"] == []
    linear = discover_linear_observables_from_events(model, ["E_LR"])
    assert linear["status"] == "CLOSED"
    assert linear["basis"] == ["E_LR"]
    assert linear["pending"] == linear["obligations"] == []


def test_zero_generator_needs_no_expansion_budget():
    from mfobs.core.linear_closure import close_linear_span

    def unexpected_image(*args):
        raise AssertionError('The zero generator has no channel images')

    result = close_linear_span(seeds={'X': {'X': 1}}, channels=(),
                               image=unexpected_image, maximum_expansions=0)
    assert result.status == 'CLOSED'
    assert result.pending == result.obligations == []


def test_zero_event_numerics_keep_the_initial_counts(invariant_query):
    from mfobs.numerics.events import affine_generator, integrate_events, simulate_events

    backend, inputs = invariant_query
    record = backend.run(inputs, uniform_matches=True)
    names, matrix, offset = affine_generator(record, {})
    assert names == ("E_LR",)
    assert np.array_equal(matrix, [[0]]) and np.array_equal(offset, [0])
    times = np.linspace(0, 1, 3)
    samples, trace = simulate_events(record, {"E_LR": 3}, {}, times, np.random.default_rng(1))
    assert np.array_equal(samples, [[3], [3], [3]])
    assert trace.empty and trace.attrs["event_count"] == 0
    moments = integrate_events(record, {"E_LR": 3}, {}, times, query="E_LR")
    assert np.allclose(moments.mean, 3) and np.allclose(moments.variance, 0)


def test_zero_event_presentations_are_closed_not_unknown(invariant_query):
    from mfobs.formatting import render
    from notebooks.paper_view import show_events

    backend, inputs = invariant_query
    record = backend.run(inputs, uniform_matches=True)
    text = render(record)
    assert "Event classes: 0" in text and "D = 0" in text
    assert "No finite event table" not in text
    assert show_events(record).empty
