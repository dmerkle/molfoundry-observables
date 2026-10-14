"""Direct molecular operands, independent chemistry and explicit failure boundaries."""
from mfobs.core.graph import edge_rows, make_graph, node_labels, node_rows
from mfobs import analyze, MolecularROP

import mfobs.models.molecular_rop.events as molecular_events
import mfobs.models.molecular_rop.inputs as molecular_inputs
from dataclasses import replace
import pytest
import sympy as sp

from mfobs.core.engine import mf

from mfobs.core.gml import disjoint_union
from mfobs.core.pattern_support import UnsupportedReduction
from mfobs.models.molecular_rop import counting as molecular
from mfobs.checks import molecular_rop as checks
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.aba import counting as terminal


@pytest.fixture(scope="module")
def baseline():
    return molecular.derive("baseline",("B",))


def test_direct_composition_does_not_call_typed_source_or_query_constructors(monkeypatch):
    def forbidden(*args,**kwargs):
        pytest.fail("typed constructor used during molecular discovery")
    monkeypatch.setattr(terminal,"source_rules",forbidden)
    monkeypatch.setattr(terminal,"query_pattern",forbidden)
    selection=rop.configured_variant(motifs=("ABA",),gated=False)
    support,sources=molecular.derive(selection,("ABA",))
    assert support.status=="POLYNOMIAL_SUPPORT_CLOSED"
    assert set(support.patterns)=={"C_ABA","M_A","M_B","S_AB","S_A"}
    assert support.images["C_ABA"]["insert_A"]==sp.Symbol("M_A")*sp.Symbol("S_AB")
    trace=[t for t in support.trace if t["kind"]=="generator-image"]
    assert trace and all(t["include_empty"] for t in trace)
    for entry in trace:
        for gml in entry["actual_operands"].values():
            assert 'label "C"' in gml
            assert 'label "Position"' not in gml
    assert any(t.get("terms") for t in trace)


def test_molecular_rule_application_counts_and_full_joint_jumps(baseline):
    support,sources=baseline
    model=molecular_events.polynomial_jump_model(support,sources,name="test",input_definition_sha256="test")
    validation=checks.validate_molecular_hosts("baseline",support,sources,jump_model=model,maximum_length=4)
    assert validation["status"]=="PASS"
    assert validation["states"]==15
    assert validation["individual_molecular_successors"]==30
    comparison=checks.compare_typed_reference(support,sources,("B",))
    assert comparison["status"]=="PASS"
    assert comparison["joint_jump_rate_identities"]==2


def test_changed_molecular_source_cannot_borrow_uniformity_certificate(baseline):
    support,sources=baseline
    changed=(replace(sources[0],rate="changed"),*sources[1:])
    assert molecular_events.compile_uniform_events(support,changed)["status"]=="UNKNOWN_UNIFORM_MATCH_OBLIGATION"


def test_zero_context_certificate_and_unknown_positive_component():
    reducer=molecular.MolecularReducer("baseline")
    triangle=make_graph({0:"C",1:"C",2:"C"},[(0,1,"-"),(1,2,"-"),(0,2,"-")])
    zero=reducer.reduce(molecular_inputs._identity("triangle",triangle))
    assert zero.expression==0
    assert zero.evidence["reduction"]["kind"]=="no-positive-embedding-in-complete-local-context-cover"
    carbonyl=make_graph({0:"C",1:"O"},[(0,1,"=")])
    with pytest.raises(UnsupportedReduction,match="positive molecular component"):
        reducer.reduce(molecular_inputs._identity("partial",carbonyl))


def test_context_cache_reuses_hosts_without_retaining_finished_reducers():
    import gc
    import weakref

    reducer = molecular.MolecularReducer("baseline")
    contexts = reducer._context_hosts(2)
    assert reducer._context_hosts(2) is contexts
    reference = weakref.ref(reducer)
    del reducer
    gc.collect()
    assert reference() is None


def test_injective_component_products_are_not_assumed_independent():
    reducer=molecular.MolecularReducer("baseline")
    two_pools=disjoint_union([rop.lactone_core("B")]*2)
    pool = sp.Symbol("M_B")
    assert reducer.reduce(molecular_inputs._identity("two",two_pools)).expression == sp.expand(pool*(pool-1))
    marker=disjoint_union([rop.reservoir("B")]*2)
    assert reducer.reduce(molecular_inputs._identity("duplicate-marker",marker)).expression==0


def test_molecular_gate_absence_and_contradiction():
    reducer=molecular.MolecularReducer("noedge")
    open_gate=molecular_inputs._identity("open",rop.catalyst(False),[(0,4)])
    assert reducer.reduce(open_gate).expression==sp.Symbol("X")
    closed_gate=molecular_inputs._identity("closed",rop.catalyst(True))
    assert reducer.reduce(closed_gate).expression==1-sp.Symbol("X")
    contradiction=molecular_inputs._identity("impossible",rop.catalyst(True),[(0,4)])
    assert reducer.reduce(contradiction).expression==0


def test_deleted_vertices_and_unhandled_conditions_are_refused():
    reducer=molecular.MolecularReducer("baseline")
    deleted=mf.Rule.fromGMLString('rule [ left [ node [ id 0 label "C" ] ] ]',add=False)
    with pytest.raises(UnsupportedReduction,match="vertex deletion"):
        reducer.reduce(deleted)
    condition=mf.Rule.fromGMLString('rule [ context [ node [ id 0 label "C" ] ] constrainAdj [ id 0 op "=" count 1 nodeLabels [ label "O" ] edgeLabels [ label "=" ] ] ]',add=False)
    with pytest.raises(UnsupportedReduction,match="condition"):
        reducer.reduce(condition)


def test_incomplete_standard_route_preserves_pending_work():
    backend=MolecularROP('baseline')
    result=analyze(backend.specification(),backend=backend,support_expansion_limit=0).record()
    assert result["status"]=="UNKNOWN_EXPANSION_LIMIT"
    assert result["support_discovery"]["pending"]
    assert result["pjc"]["status"]=="UNKNOWN_PREREQUISITE"
    assert "molecular_graph_validation" not in result
    assert "operator" not in result["pjc"]


def test_no_typed_projection_is_claimed_in_standard_contract():
    backend=MolecularROP('baseline')
    result=analyze(backend.specification(),backend=backend,linear_limit=24).record()
    assert result["status"]=="ANALYZED"
    assert result["linear"]["status"]=="CLOSED"
    assert len(result["linear"]["basis"])==3
    assert result["pjc"]["event_class_count"]==2
    from mfobs.checks.molecular_rop import validate_result
    checked = validate_result(result)
    assert checked["molecular_graph_validation"]["status"]=="PASS"
    assert checked["typed_reference_comparison"]["joint_jump_rate_identities"]==2
    assert result["composition"]["projection"]["kind"]=="no typed projection before composition"


def test_molecular_composition_budget_retains_exact_operands_and_pending_query():
    support,sources=molecular.derive("baseline",("B",),maximum_candidates=0)
    assert support.status=="UNKNOWN_COMPOSITION_LIMIT"
    assert "N_B" in support.pending
    failure=support.trace[-1]["minimal_reproducer"]
    assert failure["kind"]=="candidate-overlap-limit"
    assert "query_gml" in failure and "source_gml" in failure
    assert molecular_events.compile_uniform_events(support,sources)["status"]=="UNKNOWN_SUPPORT"


def test_profile_pruning_keeps_overlaps_repaired_by_future_identifications():
    from mfobs.core import engine
    bond=molecular_inputs._identity("bond",make_graph({0:"C",1:"H"},[(0,1,"-")]))
    spec={"atomLabels":["C","H"],"allowedBondPairs":{"-":[["C","H"]]},
          "bondOrders":{"-":1},"maxValence":{"C":1,"H":1},"incidentBounds":[],"strictSignature":True}
    # Gluing just C would have two H neighbours, but subsequently identifying
    # H repairs that overlap. A pruning rule must retain the complete bond map.
    expected=engine.grouped_restricted_common_compositions(bond,bond,forbidden_specification=spec)
    spec["allowedIncidentProfiles"]={"C":[{"-|H":1}],"H":[{"-|C":1}]}
    actual=engine.grouped_restricted_common_compositions(bond,bond,forbidden_specification=spec)
    assert actual[1]==expected[1]
    assert sum(actual[1].values())==2  # full overlap and empty overlap


def test_declared_radius_one_profiles_cover_longer_mixed_chains():
    selection=rop.configured_variant(motifs=("ABA",),gated=True)
    profiles=molecular_inputs.forbidden_specification(selection)["allowedIncidentProfiles"]
    from collections import Counter
    for word in ((),tuple("AB"),tuple("BA"),tuple("ABABA"),tuple("BBAABB")):
        graph=rop.chain(word);labels=node_labels(graph)
        for atom,label in node_rows(graph):
            ports=dict(Counter(bond+"|"+labels[v if u==atom else u] for u,v,bond in edge_rows(graph) if atom in (u,v)))
            assert ports in profiles[label]
