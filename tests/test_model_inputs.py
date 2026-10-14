"""The explicit inputs must denote exactly the analyzed quantities."""
from dataclasses import replace
from fractions import Fraction

import molfoundry as mf
import pytest

from mfobs import (BrusselatorROP, Chlorination, DNAWalker, EightSlots, GrowingRing,
                   ABA, Observable, RemyTrees, Ring8, analyze,
    MolecularROP,
)
from mfobs.state_space import StateSpace


@pytest.mark.parametrize("application", [
    ABA.motif('ABA'), ABA.motif('ABA', gated=True),
    MolecularROP('baseline'),
    Chlorination(), Chlorination("ablated"), RemyTrees(), EightSlots(),
    Ring8(), GrowingRing(), BrusselatorROP(), DNAWalker(),
])
def test_declarations_are_native_and_contain_only_selected_queries(application):
    inputs = application.specification()
    assert inputs.omega.kind and inputs.omega.definition
    assert inputs.rules and inputs.observables
    for declaration in inputs.rules.values():
        assert isinstance(declaration.rule, mf.Rule)
        inputs.omega.validate_rule(declaration.rule)
    for name, query in inputs.observables.items():
        assert name == query.name
        for _, identity in query.identity_terms():
            assert isinstance(identity, mf.Rule)
            inputs.omega.validate_rule(identity)
    if isinstance(application, ABA) and application.motifs == ("ABA",):
        assert set(inputs.observables) == {"C_ABA"}


def test_remy_query_is_one_sum_of_three_identities():
    inputs = RemyTrees(3).specification()
    assert list(inputs.observables) == ["P3"]
    assert len(inputs.observables["P3"].identity_terms()) == 3


@pytest.mark.parametrize("field", ["omega", "rule", "rate", "factor", "identifier", "observable"])
def test_changed_inputs_are_rejected_before_backend_execution(field, monkeypatch):
    app = Chlorination()
    inputs = app.specification()
    if field == "omega":
        changed = replace(inputs, omega=StateSpace(inputs.omega.kind, {**inputs.omega.definition, "chlorine_total": 5}, inputs.omega.chemistry))
    elif field == "observable":
        query = inputs.observables["Z"]
        _, identity = query.identity_terms()[0]
        changed = replace(inputs, observables={"Z": Observable("Z", terms=((Fraction(2), identity),))})
    else:
        rules = dict(inputs.rules)
        declaration = rules["chlorinate"]
        kwargs = ({"rate": "another_rate"} if field == "rate" else
                  {"per_match": Fraction(2)} if field == "factor" else
                  {"identifier": "another_rule"} if field == "identifier" else
                  {"rule": inputs.rules["grow"].rule})
        rules["chlorinate"] = replace(declaration, **kwargs)
        changed = replace(inputs, rules=rules)
    def forbidden_run(*args, **kwargs):
        pytest.fail("unsupported input reached analysis")
    from mfobs.models.chlorination.counting import LocalConditionalReducer
    monkeypatch.setattr(LocalConditionalReducer, "activity", forbidden_run)
    with pytest.raises(ValueError, match="unsupported Ω"):
        analyze(changed, backend=app)


def test_actual_omega_is_forwarded_to_backend(monkeypatch):
    from mfobs.models.remy import analysis as remy_trees
    inputs = RemyTrees().specification()
    def run(**kwargs):
        assert kwargs["state_space"] is inputs.omega
        return {"linear": {"status": "CLOSED"}}
    monkeypatch.setattr(remy_trees, "common_analysis", run)
    assert analyze(inputs, backend=RemyTrees()).linear.status == "CLOSED"


def test_brusselator_count_expressions_count_molecules():
    from mfobs.core.gml import graph_gml
    from mfobs.models import brusselator_growth_chemistry as chemistry
    from mfobs.models.molecular_rop import inputs as copolymer_rop
    inputs = BrusselatorROP(("X", "Y", "M", "N")).specification()
    hosts = {"X": chemistry.species("X"), "Y": chemistry.species("Y"),
             "M": copolymer_rop.monomer("B"), "N": copolymer_rop.chain(("B", "B"))}
    for name, host in hosts.items():
        coefficient, identity = inputs.observables[name].identity_terms()[0]
        count = identity.numMatches([mf.Graph.fromGMLString(graph_gml(host), add=False)], mf.ApplyMode.UseAny)
        assert coefficient * count == (2 if name == "N" else 1)
    assert inputs.rules["autocatalysis"].per_match == Fraction(1, 27648)


def test_guarded_and_ablated_share_omega_and_query():
    guarded, ablated = Chlorination().specification(), Chlorination("ablated").specification()
    assert guarded.omega.record() == ablated.omega.record()
    assert guarded.signature()[2] == ablated.signature()[2]
    assert guarded.rules["grow"].rule.canonKey == ablated.rules["grow"].rule.canonKey
    assert guarded.rules["chlorinate"].rule.matchCondition is not None
    assert ablated.rules["chlorinate"].rule.matchCondition is None
