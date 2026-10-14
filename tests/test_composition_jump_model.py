"""Model-independent checks for composition-derived jump compilation."""
import inspect

from mfobs.core.composition_jump_model import (
    compare_polynomial_jump_models,
    compile_composition_jump_model,
)
from mfobs.core.recursive_words import TerminalGrowthGrammar
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.core.typed_incidence import dangling_certificate, typed_incidence_commutator


def test_compiler_api_cannot_accept_a_supplied_delta_table():
    parameters = inspect.signature(compile_composition_jump_model).parameters
    assert "delta" not in parameters
    assert "events" not in parameters
    assert {
        "observables", "reactions", "pair_commutator", "activity_reducer",
        "jump_closure", "jump_projection_name",
    } <= set(parameters)


def test_nonchemical_alphabet_uses_the_same_composition_pipeline():
    grammar = TerminalGrowthGrammar(
        "generic-two-letter-word",
        ("u", "v"),
        (("u", "v"),),
        True,
        "0" * 64,
    )
    graph_model = TerminalWordGraphModel(grammar)
    specification = graph_model.specification()

    def pair(observable, reaction):
        return typed_incidence_commutator(
            observable,
            reaction,
            specification=specification,
            include_empty=False,
        )

    compilation = compile_composition_jump_model(
        name=grammar.name,
        ring=grammar.polynomial_ring(),
        observables=graph_model.observable_expressions(),
        reactions=graph_model.reaction_classes(),
        pair_commutator=pair,
        activity_reducer=graph_model.activity_polynomial,
        event_partition_certificate=graph_model.event_partition_certificate(),
        input_definition_sha256=grammar.input_definition_sha256,
        dangling_certificate=dangling_certificate(specification),
    )
    comparison = compare_polynomial_jump_models(
        compilation.model,
        grammar.polynomial_jump_model(),
    )
    assert comparison["status"] == "PASS"
    assert len(compilation.model.events) == 10
    assert compilation.evidence["backend"] == (
        "MolFoundry DPO common-overlap composition"
    )
    assert all(
        event["observable_commutators"]
        for event in compilation.evidence["event_classes"]
    )
    assert compilation.model.provenance.record()["supplied_delta_table"] is False
