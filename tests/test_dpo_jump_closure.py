"""DPO jump closure I <- K -> I (Behr et al., arXiv:2106.02573, Thm. 6(2), Remark 3)."""
from mfobs.core.graph import edge_rows, node_labels, make_graph
from mfobs.core.rules import make_rule
import pytest

from mfobs import Observable, RatedRule, SimpleGraphs, SimpleGraphBackend, ModelInputs, analyze


from mfobs.core.engine import dpo_jump_closure, mf, sqpo_jump_closure
from mfobs.core.gml import rule_side
from mfobs.core.identity import DanglingNotCertified, certified_input_identity, jump_input




def applications(rule, host):
    """DPO applications of ``rule`` to ``host``: composites of (create host; rule) with empty input."""
    create = make_rule("host", make_graph({}), host)
    expression = create * mf.rcCommon(connected=False, maximum=False, includeEmpty=True) * rule
    composites = mf.RCEvaluator([create, rule]).eval(expression, onlyUnique=False)
    return sum(1 for c in composites if not node_labels(rule_side(c, "left")))


DELETE_A = make_rule("delete", make_graph({0: "A"}), make_graph({}))
ISOLATED = make_graph({0: "A"})
ATTACHED = make_graph({0: "A", 1: "B"}, [(0, 1, "-")])


def test_dpo_jump_closure_counts_dpo_applications_including_dangling():
    jump = dpo_jump_closure(DELETE_A)

    assert [applications(DELETE_A, h) for h in (ISOLATED, ATTACHED)] == [1, 0]
    assert [applications(jump, h) for h in (ISOLATED, ATTACHED)] == [1, 0]
    # The SqPO closure I <- I -> I has no dangling condition.
    assert [applications(sqpo_jump_closure(DELETE_A), h) for h in (ISOLATED, ATTACHED)] == [1, 1]
    assert [v.membership.value for v in jump.vertices] == ["deleted", "created"]


def test_jump_input_reads_deleted_vertices_and_rejects_other_spans():
    parsed, deleted = jump_input(dpo_jump_closure(DELETE_A))
    assert (dict(parsed.labels), deleted) == ({0: "A"}, (0,))
    assert jump_input(make_rule("id", ATTACHED, ATTACHED))[1] == ()
    with pytest.raises(ValueError):
        jump_input(DELETE_A)


def test_public_condition_api_preserves_nested_scope_and_fresh_ids():
    from molfoundry.constraints import ConditionNot, ConditionExists
    rule = mf.Rule.fromGMLString('''rule [
      left [ node [ id 10 label "A" ] ]
      condition [ not [ exists [ node [ id 11 label "X" ]
        not [ exists [ node [ id 12 label "Y" ]
          edge [ source 11 target 12 label "-" ] ] ] ] ] ]
    ]''', add=False)
    assert isinstance(rule.matchCondition, ConditionNot)
    assert isinstance(rule.matchCondition.condition, ConditionExists)
    jump = dpo_jump_closure(rule)
    assert jump.matchCondition == rule.matchCondition
    assert {v.id for v in jump.vertices if v.isCreated}.isdisjoint({10, 11, 12})
    for labels, count in [({0: 'A'}, 1), ({0: 'A', 1: 'X'}, 0)]:
        from mfobs.core.gml import graph_gml
        host = mf.Graph.fromGMLString(graph_gml(make_graph(labels)), add=False)
        for candidate in (rule, jump):
            dg = mf.DG(graphDatabase=[host])
            with dg.build() as builder:
                assert len(builder.apply([host], candidate, mode=mf.ApplyMode.UseAny, verbosity=0)) == count


def test_public_span_preserves_context_edges_and_relabelling():
    from mfobs.core.rules import rule_graph, rule_vertex_ids
    from mfobs.core.graph import same_graph
    rule = mf.Rule.fromGMLString('''rule [
      context [ node [ id 10 label "A" ] node [ id 20 label "B" ]
        edge [ source 10 target 20 label "-" ] ]
      left [ node [ id 30 label "X" ] ]
      right [ node [ id 30 label "Y" ] ]
    ]''', add=False)
    assert edge_rows(rule_graph(rule, "context")) == ((10, 20, '-'),)
    assert rule_vertex_ids(rule, "changed") == (30,)
    jump = dpo_jump_closure(rule)
    assert same_graph(rule_graph(jump), rule_graph(jump, "right"))
    assert same_graph(rule_graph(jump), rule_graph(rule))


def test_input_identity_needs_a_dangling_certificate():
    jump = dpo_jump_closure(DELETE_A)
    with pytest.raises(DanglingNotCertified):
        certified_input_identity(jump)

    def refuse(graph, deleted):
        raise DanglingNotCertified("no")

    with pytest.raises(DanglingNotCertified):
        certified_input_identity(jump, refuse)
    reduced = certified_input_identity(jump, lambda graph, deleted: None)
    assert reduced.canonKey == make_rule("id", ISOLATED, ISOLATED).canonKey


def test_typed_incidence_certificate_requires_saturated_ports():
    from mfobs.core.typed_incidence import dangling_certificate

    certify = dangling_certificate({
        "objectLabels": ["V"], "edgeNodeLabels": ["E"],
        "sourceBond": "src", "targetBond": "tgt", "atMostOneOutgoing": False,
    })
    complete = make_graph({0: "V", 1: "E", 2: "V"}, [(0, 1, "src"), (1, 2, "tgt")])
    certify(complete, (1,))
    with pytest.raises(DanglingNotCertified):
        certify(make_graph({0: "V", 1: "E"}, [(0, 1, "src")]), (1,))
    with pytest.raises(DanglingNotCertified):
        certify(complete, (0,))  # an object may have further incidences


def test_simple_graphs_do_not_drop_the_dangling_condition():
    omega = SimpleGraphBackend(SimpleGraphs(vertex_labels={"L", "R"}, edge_labels={"-"}))
    delete = RatedRule(mf.ruleGMLString('''rule [
      ruleID "delete"
      context [ node [ id 0 label "L" ] node [ id 2 label "L" ] ]
      left [ node [ id 1 label "R" ]
        edge [ source 0 target 1 label "-" ]
        edge [ source 1 target 2 label "-" ] ]
    ]''', add=False), "k3")
    query = Observable("N_R", mf.graphGMLString(
        'graph [ node [ id 0 label "R" ] ]', add=False))
    result = analyze(ModelInputs((omega).omega, {r.identifier:r for r in [delete]}, {(query).name:query}), backend=omega, uniform_matches=True)

    assert result.support.status == "UNKNOWN_UNSUPPORTED_REDUCTION"
    assert not result.record()["events"]
    inputs = ModelInputs(omega.omega, {delete.identifier: delete}, {query.name: query})
    native_linear = omega.linear_closure(inputs)
    step = native_linear.steps[0]
    (name, coefficient), = step["coordinates"].items()
    assert coefficient == "-1"
    # On a star the third edge dangles: delete never applies, so G N_R = 0 there.
    star = make_graph({0: "R", 1: "L", 2: "L", 3: "L"}, [(0, 1, "-"), (0, 2, "-"), (0, 3, "-")])
    observable = native_linear.basis[name]
    assert sum(c * applications(observable.representatives[k], star)
               for k, c in observable.coefficients.items()) == 0
