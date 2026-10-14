"""Native MolFoundry inputs for the graph tutorial and SimpleGraphs backend."""

from pathlib import Path
from runpy import run_path

import pytest
from molfoundry import graphGMLString, ruleGMLString

from mfobs import Observable, RatedRule, SimpleGraphs, SimpleGraphBackend, ModelInputs, analyze


def test_public_simple_graphs_tutorial_runs_both_algorithms():
    tutorial = Path(__file__).parents[1] / "tutorials/graph_closure_tutorial.py"
    namespace = run_path(str(tutorial))
    linear, jump = namespace["linear"], namespace["jump"]

    assert linear.status == "CLOSED"
    assert list(linear.basis) == ["N_L", "E_LR"]
    assert list(jump.support.patterns) == ["N_L", "E_LR"]

    event, = jump.compilation.model.events
    assert event.channel == "k"
    assert event.delta == (1, 0)
    assert jump.compilation.model.ring.text(event.propensity) == "E_LR"
    assert jump.pjc["status"] == "CLOSED"


def test_native_context_edge_and_nonconsecutive_ids():
    rule = ruleGMLString('''rule [
      ruleID "birth_at_edge"
      context [
        node [ id 7 label "L" ] node [ id 13 label "R" ]
        edge [ source 7 target 13 label "-" ]
      ]
      right [ node [ id 29 label "L" ] ]
    ]''', add=False)
    query = graphGMLString('graph [ node [ id 42 label "L" ] ]', add=False)
    omega = SimpleGraphs(vertex_labels={"L", "R"}, edge_labels={"-"})
    inputs = ModelInputs(omega, {"birth_at_edge": RatedRule(rule, rate="b")},
                         {"N_L": Observable("N_L", query)})
    result = analyze(inputs, backend=SimpleGraphBackend(omega), uniform_matches=True)

    assert result.linear.status == result.pjc.status == "CLOSED"
    assert list(result.support.data["patterns"]) == ["N_L", "E_LR"]
    event, = result.record()["events"]
    assert event["rate"] == "b"
    assert event["jump"] == (1, 0)
    assert event["activity"]["text"] == "E_LR"


@pytest.mark.parametrize("identifiers", [(None,), ("alias_birth",), ("first", "second")])
def test_simple_graphs_preserves_declared_rule_identifiers(identifiers):
    rule = ruleGMLString('''rule [ ruleID "native_birth"
        context [node [id 0 label "L"]]
        right [node [id 1 label "L"]] ]''', add=False)
    omega = SimpleGraphs(vertex_labels={"L"}, edge_labels={"-"})
    rules = {identifier or rule.name: RatedRule(rule, rate=f"k_{index}", identifier=identifier)
             for index, identifier in enumerate(identifiers)}
    query = Observable("N_L", graphGMLString('graph[node[id 0 label "L"]]', add=False))
    inputs = ModelInputs(omega, rules, {query.name: query})
    result = SimpleGraphBackend(omega).jump_closure(inputs, uniform_matches=True)

    assert result.support.status == "POLYNOMIAL_SUPPORT_CLOSED"
    assert list(result.support.activities) == list(rules)
    assert list(result.support.images["N_L"]) == list(rules)
    assert [event["source_rule"] for event in result.compilation.evidence["events"]] == list(rules)
    assert [item["source_rule"] for item in result.compilation.evidence["uniform_match_obligations"]] == list(rules)
    for index, (identifier, event) in enumerate(zip(rules, result.compilation.model.events)):
        assert event.identifier == f"{identifier}-support-{index}"
        assert event.rule == identifier
        assert event.channel == f"k_{index}"
        assert event.delta == (1,)
        assert result.compilation.model.ring.text(event.propensity) == "N_L"
    assert rule.name == "native_birth"


def test_simple_graph_input_digest_tracks_aliases_and_ignores_their_order():
    rule = ruleGMLString('''rule [ ruleID "native_birth"
        context [node [id 0 label "L"]]
        right [node [id 1 label "L"]] ]''', add=False)
    omega = SimpleGraphs(vertex_labels={"L"}, edge_labels={"-"})
    query = Observable("N_L", graphGMLString('graph[node[id 0 label "L"]]', add=False))
    backend = SimpleGraphBackend(omega)

    def digest(identifiers):
        rules = {name: RatedRule(rule, identifier=name) for name in identifiers}
        result = backend.jump_closure(ModelInputs(omega, rules, {query.name: query}),
                                      uniform_matches=True)
        return result.compilation.model.input_definition_sha256

    assert digest(("first",)) != digest(("second",))
    assert digest(("first", "second")) == digest(("second", "first"))


@pytest.mark.parametrize("condition", [
    'constrainNoEdge [ source 0 target 1 ]',
    'constrainAdj [ id 0 op "=" count 0 nodeLabels [ label "R" ] ]',
    '''condition [ not [ exists [ node [ id 7 label "R" ]
         edge [ source 0 target 7 label "-" ]
         not [ exists [ node [ id 8 label "L" ]
           edge [ source 7 target 8 label "-" ] ] ] ] ] ]''',
])
def test_simple_graphs_rejects_conditions_instead_of_dropping_them(condition):
    rule = ruleGMLString('''rule [
      context [ node [ id 0 label "L" ] node [ id 1 label "R" ] ]
      right [ node [ id 2 label "L" ] ]
    ''' + condition + ']', add=False)
    query = graphGMLString('graph [ node [ id 0 label "L" ] ]', add=False)
    omega = SimpleGraphs(vertex_labels={"L", "R"}, edge_labels={"-"})
    inputs = ModelInputs(omega, {"birth": RatedRule(rule)},
                         {"N_L": Observable("N_L", query)})
    with pytest.raises(ValueError, match="unconditional rules only"):
        analyze(inputs, backend=SimpleGraphBackend(omega))


def test_simple_graph_input_digest_is_independent_of_backend_instance():
    rule = ruleGMLString('''rule [
      ruleID "birth"
      context [ node [ id 0 label "L" ] ]
      right [ node [ id 1 label "L" ] ]
    ]''', add=False)
    query = Observable("N_L", graphGMLString(
        'graph [ node [ id 0 label "L" ] ]', add=False))
    omega = SimpleGraphs(vertex_labels={"L", "R"}, edge_labels={"-"})
    inputs = ModelInputs(omega, {"birth": RatedRule(rule)}, {query.name: query})
    first, second = SimpleGraphBackend(omega), SimpleGraphBackend(omega)

    def digest(backend, model_inputs):
        result = backend.jump_closure(model_inputs, uniform_matches=True)
        return result.compilation.model.input_definition_sha256

    assert digest(first, inputs) == digest(second, inputs)
    larger_omega = SimpleGraphs(vertex_labels={"L", "R", "S"}, edge_labels={"-"})
    larger_inputs = ModelInputs(larger_omega, inputs.rules, inputs.observables)
    assert digest(first, inputs) != digest(SimpleGraphBackend(larger_omega), larger_inputs)


def test_simple_graph_input_digest_is_independent_of_rule_mapping_order():
    omega = SimpleGraphs(vertex_labels={"L"}, edge_labels={"-"})
    rules = {name: RatedRule(ruleGMLString(f'''rule [ ruleID "{name}"
        context [node [id 0 label "L"]]
        right [node [id 1 label "L"]] ]''', add=False), rate=name)
        for name in ("first", "second")}
    query = Observable("N_L", graphGMLString('graph[node[id 0 label "L"]]', add=False))
    backend = SimpleGraphBackend(omega)
    inputs = [ModelInputs(omega, ordering, {query.name: query})
              for ordering in (rules, dict(reversed(tuple(rules.items()))))]
    assert inputs[0].signature() == inputs[1].signature()
    results = [backend.jump_closure(model, uniform_matches=True) for model in inputs]
    assert results[0].compilation.model.input_definition_sha256 == results[1].compilation.model.input_definition_sha256
