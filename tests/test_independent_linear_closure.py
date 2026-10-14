"""Algorithm 1 works without Algorithm 2 and preserves joint jump effects."""
from unittest.mock import patch

import pytest
import sympy as sp

from mfobs import BrusselatorROP, Observable, ABA, MolecularROP
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.core.gml import disjoint_union
from mfobs.core.polynomial_jump_closure import generator_action
from mfobs.models.aba import counting, events
from mfobs.models.molecular_rop import inputs as molecular
from mfobs.models.aba import analysis as path_analysis
from mfobs.models.molecular_rop import analysis as molecular_analysis
from mfobs.models import brusselator_pjc


def test_ambiguous_source_identifiers_are_rejected_before_composition():
    rules = counting.source_rules()
    with pytest.raises(ValueError, match="identifiers must be unique"):
        discover_linear_observables((rules[0], rules[0]),
            {"C_ABA": counting.query_pattern("ABA")}, counter=counting.PathCounter())


@pytest.mark.parametrize("representation", ["path", "molecular"])
def test_linear_api_needs_neither_support_nor_events(representation):
    backend = (MolecularROP if representation == "molecular" else ABA).motif('ABA')
    with (patch("mfobs.models.aba.counting.discover_support", side_effect=AssertionError("Algorithm 2")),
          patch("mfobs.models.molecular_rop.counting.discover_support", side_effect=AssertionError("Algorithm 2")),
          patch("mfobs.models.aba.events.polynomial_jump_model", side_effect=AssertionError("events"))):
        result = backend.linear_closure(backend.specification(), maximum_dimension=4)
    assert result["status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert result["basis"][0] == "C_ABA"
    assert result["steps"][2]["generator_image"]["text"] == "M_A*S_AB"


@pytest.mark.parametrize("representation", ["path", "molecular"])
def test_jump_api_does_not_run_algorithm_one(representation):
    backend = (MolecularROP if representation == "molecular" else ABA).motif('ABA')
    with patch.object(path_analysis, "discover_linear_observables",
               side_effect=AssertionError("Algorithm 1")), patch.object(
               molecular_analysis, "discover_linear_observables", side_effect=AssertionError("Algorithm 1")):
        result = backend.jump_closure(backend.specification())
    assert result["linear"] is None
    assert result["pjc"]["status"] == "CLOSED"


@pytest.mark.parametrize("representation", ["path", "molecular"])
def test_jump_search_is_not_limited_by_the_linear_basis_budget(representation):
    selection = molecular.configured_variant(motifs=('ABA','ABBA'), gated=False)
    backend = (MolecularROP(selection) if representation == "molecular" else
               ABA(tuple(''.join(word) for word in selection.motifs),
                   alphabet=selection.alphabet, gated=selection.gated))
    result = backend.jump_closure(backend.specification(), linear_limit=1)
    assert result['linear'] is None
    assert result['pjc']['status'] == 'CLOSED'


def test_brusselator_linear_api_does_not_compile_events():
    backend = BrusselatorROP()
    with patch.object(brusselator_pjc, "composition_compilation", side_effect=AssertionError("events")):
        linear = backend.linear_closure(backend.specification(), maximum_dimension=6)
    assert linear["status"] == "UNKNOWN_DIMENSION_LIMIT"


@pytest.mark.parametrize("queries, limit", [(('X','Y'),24), (('M',),6), (('N',),6)])
def test_brusselator_direct_images_match_the_separate_event_generator(queries, limit):
    linear = brusselator_pjc.linear_closure(queries, linear_limit=limit)
    assert linear["status"] == "UNKNOWN_DIMENSION_LIMIT"
    model = brusselator_pjc.polynomial_jump_model()
    for step in linear["steps"]:
        source = model.ring.add(*[model.ring.scale(model.ring.monomial(t["monomial"]), t["coefficient"])
            for t in linear["basis_expressions"][step["source"]]["terms"]])
        assert decode(step["generator_image"]) == decode(model.ring.record(
            generator_action(model, source, channel=step["channel"])))


def decode(record):
    return sp.expand(sum(sp.Rational(term["coefficient"]) *
                         sp.prod(sp.Symbol(name)**power for name, power in term["monomial"].items())
                         for term in record["terms"]))


@pytest.mark.parametrize("gated", [False, True])
@pytest.mark.parametrize("representation", ["path", "molecular"])
def test_all_direct_images_agree_with_separately_derived_event_generator(gated, representation):
    backend = (MolecularROP if representation == 'molecular' else ABA).motif('ABA', gated=gated)
    linear = backend.linear_closure(backend.specification(), maximum_dimension=24)
    assert linear["status"] == "UNKNOWN_DIMENSION_LIMIT"
    if representation == "molecular":
        from mfobs.models.molecular_rop import counting as molecule_counting, events as molecule_events
        support, sources = molecule_counting.derive(backend.configuration(), ("ABA",))
        compiler = molecule_events.polynomial_jump_model
    else:
        support, sources = counting.derive(("ABA",), gated=gated)
        compiler = events.polynomial_jump_model
    model = compiler(support, sources, name="reference")
    ring = model.ring
    for step in linear["steps"]:
        expression = linear["basis_expressions"][step["source"]]
        source = ring.add(*[ring.scale(ring.monomial(t["monomial"]), t["coefficient"])
                            for t in expression["terms"]])
        actual = generator_action(model, source, channel=step["channel"])
        assert decode(step["generator_image"]) == decode(ring.record(actual))


def test_second_moment_is_not_computed_by_the_drift_product_rule():
    # Two disjoint B molecules count M_B(M_B-1); adding M_B gives M_B².
    one = molecular._identity("one", molecular.lactone_core("B"))
    two = molecular._identity("two", disjoint_union([molecular.lactone_core("B")]*2))
    import molfoundry as mf
    from mfobs.core.gml import graph_gml
    host = mf.Graph.fromGMLString(graph_gml(disjoint_union([molecular.monomer('B')]*3)), add=False)
    assert two.numMatches([host], onlyProper=True) == 6
    query = Observable("second_moment", terms=((1, two), (1, one)))
    from mfobs.models.molecular_rop.counting import MolecularReducer
    linear = discover_linear_observables(molecular.source_rules("baseline"),
        {query.name: query}, counter=MolecularReducer("baseline"), maximum_dimension=3)
    rows = {row["channel"]: decode(row["generator_image"])
            for row in linear["steps"] if row["source"] == query.name}
    M = sp.Symbol("M_B")
    assert rows["beta_B"] == 2*M + 1
    assert rows["k_B"] == -2*M**2 + M
