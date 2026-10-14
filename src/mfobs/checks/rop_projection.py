"""Suffix-refined ROP projection for independent checks and recorded comparisons."""
from __future__ import annotations
from functools import lru_cache
from mfobs.core.analysis_contract import (
    AnalysisContract,
    NormalizedRule,
    ObservableDeclaration,
    StateLanguageContract,
)
from mfobs.core.composition_jump_model import compile_composition_jump_model
from mfobs.core.gml import rule_gml
from mfobs.core.conditions import engine_no_edges
from mfobs.core.recursive_words import TerminalGrowthGrammar
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.core.typed_incidence import dangling_certificate, pair_commutator
from mfobs.models.molecular_rop.inputs import RopConfiguration, configuration, observable_rules, rules


@lru_cache(None)
def _contract_and_grammar(
    selected: RopConfiguration,
) -> tuple[AnalysisContract, TerminalGrowthGrammar]:
    provisional = TerminalGrowthGrammar(
        selected.identifier, selected.alphabet, selected.motifs,
        selected.gated, "pending",
    )
    graph_observables = observable_rules(provisional)
    normalized_rules = []
    for identifier, rule in rules(selected).items():
        rate = (
            f"beta_{identifier[-1]}" if identifier.startswith("feed_") else
            f"k_{identifier[-1]}" if identifier.startswith("insert_") else
            "a" if identifier == "donor_bind" else "b"
        )
        conditions = (
            ("NoEdge(Pd,N)",)
            if engine_no_edges(rule) else ()
        )
        normalized_rules.append(NormalizedRule(
            identifier, rate, "conditional-DPO", rule_gml(rule), conditions,
            boundary=identifier.startswith("feed_"),
        ))
    declarations = tuple(
        ObservableDeclaration(
            observable.name,
            observable.meaning,
            observable.role,
            (
                "selected scientific query and linear-discovery seed"
                if observable.role == "selected"
                else
                "state-language support coordinate; not a linear-discovery seed"
                if observable.role == "support"
                else
                "generated mechanically by the terminal-word suffix automaton"
            ),
            rule_gml(graph_observables[observable.name]),
        )
        for observable in provisional.observables()
    )
    queries = tuple(
        observable.name
        for observable in provisional.primary_observables()
        if observable.role == "selected"
    )
    language = StateLanguageContract(
        f"lactone-terminal-word-{selected.identifier}",
        (
            "one benzyl-alcohol-initiated linear polyester word over "
            f"{selected.alphabet}, arbitrary free monomer pools, persistent feed "
            + ("reservoirs, and one open/closed hemilabile catalyst" if selected.gated else "reservoirs")
        ),
        True,
        (
            "exactly one unbranched chain and one terminal OH",
            "A/B order is the covalent backbone order and is never quotiented by reversal",
            "every internal insertion preserves atoms; feeds are explicit open boundaries",
            "free monomer and incorporated-repeat motifs are disjoint by ring closure",
            *(
                ("exactly one Pd fragment; X=1 iff no Pd--N edge is present",)
                if selected.gated else ()
            ),
        ),
        (
            "reject overlaps that branch the unique chain or identify distinct chain ports",
            "reject cyclic/overvalent insertion overlaps outside the declared grammar",
            "retain injective monomer-copy multiplicity and distinct A/B channels",
            *(
                ("transport NoEdge(Pd,N) and reject positive-edge/NoEdge contradictions",)
                if selected.gated else ()
            ),
        ),
        (
            "Induction on rule applications: feed adds one declared disconnected ring; "
            "insertion consumes exactly one such ring and extends the unique terminal OH; "
            "binding toggles only the optional Pd--N edge."
        ),
        (
            "word, free-pool counts"
            + (", and Boolean open-gate state" if selected.gated else "")
        ),
    )
    contract = AnalysisContract(
        f"rop-family-{selected.identifier}", selected.title,
        tuple(normalized_rules), language, declarations, queries, queries,
        {rule.identifier: "1 per injective match" for rule in normalized_rules},
    )
    grammar = TerminalGrowthGrammar(
        selected.identifier, selected.alphabet, selected.motifs,
        selected.gated, contract.digest,
    )
    return contract, grammar


def analysis_contract(selection: str | RopConfiguration) -> AnalysisContract:
    return _contract_and_grammar(configuration(selection))[0]


def grammar(selection: str | RopConfiguration) -> TerminalGrowthGrammar:
    return _contract_and_grammar(configuration(selection))[1]


@lru_cache(None)
def _composition_compilation(selected: RopConfiguration):
    """Compile a full suffix refinement for independent checks and benchmarks."""

    language = grammar(selected)
    graph_model = TerminalWordGraphModel(language)
    specification = graph_model.specification()

    return compile_composition_jump_model(
        name=language.name,
        ring=language.polynomial_ring(),
        observables=graph_model.observable_expressions(),
        reactions=graph_model.reaction_classes(),
        pair_commutator=pair_commutator(specification, include_empty=False),
        activity_reducer=graph_model.activity_polynomial,
        event_partition_certificate=graph_model.event_partition_certificate(),
        input_definition_sha256=language.input_definition_sha256,
        dangling_certificate=dangling_certificate(specification),
    )


def composition_compilation(selection: str | RopConfiguration):
    """Compile a preset or custom configuration, cached by full model input."""

    return _composition_compilation(configuration(selection))
