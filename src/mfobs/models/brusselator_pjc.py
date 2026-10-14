"""Common-contract pJC/EMGF analysis of the Brusselator--ROP example.

The graph-level worklist in :mod:`mfobs.models.brusselator_moments` discovers
an unbounded factorial-moment hierarchy and therefore reports ``UNKNOWN``.
That fact is fully compatible with polynomial jump closure: the same six
literal DPO rules have constant jumps in a four-coordinate population quotient
and polynomial propensities.  Only X and Y are selected queries.  M and N are
mechanically generated support coordinates needed to express the monomer and
chain-growth rule classes.  This module records that second, finite EMGF
description without replacing or approximating the moment hierarchy.
"""
from __future__ import annotations
from mfobs.core.rules import make_rule

from functools import lru_cache
from fractions import Fraction

from mfobs.core.analysis_contract import (
    AnalysisContract,
    NormalizedRule,
    ObservableDeclaration,
    StateLanguageContract,
    assert_same_input,
)

from mfobs.core.engine import mf
from mfobs.core.composition_jump_model import (
    CompositionJumpCompilation,
    compile_composition_jump_model,
)
from mfobs.core.gml import rule_gml
from mfobs.core.polynomial_jump_closure import PolynomialJumpModel, PolynomialRing, analyze_polynomial_jump_closure
from mfobs.models import brusselator_growth_chemistry as chemistry
from mfobs.models.brusselator_growth_projection import (
    BrusselatorGrowthGraphModel,
)
from mfobs.models.molecular_rop import inputs as copolymer_rop
from mfobs.core.typed_incidence import pair_commutator


COORDINATES = ("X", "Y", "M", "N")
SELECTED_QUERIES = ("X", "Y")
SUPPORT_COORDINATES = ("M", "N")


def _identity(name: str, graph) -> str:
    return rule_gml(make_rule(name, graph, graph))


def coordinate_identity_gml(name: str) -> str:
    identities = {
        "X": _identity("X", chemistry.species("X")),
        "Y": _identity("Y", chemistry.species("Y")),
        "M": _identity("M", copolymer_rop.monomer("B")),
        "N": _identity("N", copolymer_rop.repeat_fragment(("B",))),
    }
    try:
        return identities[name]
    except KeyError as exc:
        raise ValueError(f"unknown Brusselator coordinate {name!r}") from exc


def selected_query_names(selection) -> tuple[str, ...]:
    """Validate common query input against the four graph coordinates."""

    known = {
        name: mf.Rule.fromGMLString(
            coordinate_identity_gml(name), add=False
        ).canonKey
        for name in COORDINATES
    }
    selected = []
    for query in selection.observables:
        if len(query.terms) != 1 or query.terms[0].coefficient != Fraction(1):
            raise ValueError("Brusselator accepts one unit identity per coordinate")
        key = mf.Rule.fromGMLString(query.terms[0].gml, add=False).canonKey
        if query.name not in known or key != known[query.name]:
            raise ValueError(
                f"{query.name}: GML is not the canonical X, Y, M or N identity"
            )
        selected.append(query.name)
    return tuple(selected)


@lru_cache(None)
def analysis_contract(selected_queries: tuple[str, ...] = SELECTED_QUERIES) -> AnalysisContract:
    graph_rules = chemistry.rules()
    boundary = {chemistry.SOURCE, chemistry.CONVERSION, chemistry.SINK,
                chemistry.MONOMER_INLET}
    rules = tuple(
        NormalizedRule(
            identifier=name,
            rate=chemistry.RATE_CHANNELS[name],
            semantics="conditional-DPO",
            gml=rule_gml(graph_rules[name]),
            boundary=name in boundary,
        )
        for name in chemistry.REACTION_ORDER
    )
    if not selected_queries or set(selected_queries) - set(COORDINATES):
        raise ValueError("Brusselator queries must be selected from X,Y,M,N")
    meanings = {
        "X": "number of 1-butene X molecules",
        "Y": "number of methylcyclopropane Y molecules",
        "M": "number of free epsilon-caprolactone molecules",
        "N": "number of incorporated caprolactone repeat units",
    }
    default_origins = {
        "X": "user-selected complete-component identity query and linear seed",
        "Y": "user-selected complete-component identity query and linear seed",
        "M": "mechanically generated from literal component inventories during composition",
        "N": "mechanically generated terminal-word coordinate during composition",
    }
    observables = tuple(
        ObservableDeclaration(
            name, meanings[name],
            "selected" if name in selected_queries else "support",
            (
                default_origins[name]
                if selected_queries == SELECTED_QUERIES else
                "user-selected graph identity and linear seed"
                if name in selected_queries else
                "mechanically available composition support coordinate"
            ),
            coordinate_identity_gml(name),
        )
        for name in COORDINATES
    )
    language = StateLanguageContract(
        identifier="brusselator-rop-populations-v1",
        carrier=(
            "all nonnegative X/Y/monomer populations and one benzyl-initiated "
            "linear polyester of arbitrary length"
        ),
        infinite=True,
        invariants=(
            "X, Y and free-monomer counts are nonnegative integers",
            "exactly one linear benzyl-initiated chain with one terminal OH",
            "N equals the number of incorporated repeat units in that chain",
            "A and B are chemostatted; D and E are open-reactor products",
        ),
        overlap_constraints=(
            "complete molecular components are pairwise disjoint",
            "two X reactants in autocatalysis are distinct injective matches",
            "growth matches one complete X, one free lactone and the unique chain end",
            "negative jumps vanish at population boundaries through their propensities",
        ),
        preservation_argument=(
            "The four core rewrites preserve molecular inventories according to their "
            "literal stoichiometry; net alcoholysis extends the unique chain by one; "
            "the declared boundary rules account explicitly for open-system material."
        ),
        semantic_quotient=(
            "population counts (X,Y,M,N), with exact ordered falling-factorial "
            "normalization for the two-X match"
        ),
    )
    return AnalysisContract(
        model_id="brusselator-growth",
        title="Bond-isomer Brusselator coupled to lactone ROP",
        rules=rules,
        state_language=language,
        observables=observables,
        linear_seeds=selected_queries,
        pjc_targets=selected_queries,
        normalization={
            "rule_matches": "automorphism-normalized injective molecular matches",
            "autocatalysis": "ordered (X)_2 Y propensity; no spurious factor 1/2",
            "parameters": "chemostat and volume factors remain in rate channels",
        },
    )


@lru_cache(None)
def composition_compilation(
    selected_queries: tuple[str, ...] = SELECTED_QUERIES,
) -> CompositionJumpCompilation:
    """Derive every production jump through generated typed DPO composition."""

    contract = analysis_contract(selected_queries)
    graph_model = BrusselatorGrowthGraphModel(
        input_definition_sha256=contract.digest
    )
    specification = graph_model.specification()

    return compile_composition_jump_model(
        name=contract.model_id,
        ring=PolynomialRing(COORDINATES),
        observables=graph_model.observable_expressions(),
        reactions=graph_model.reaction_classes(),
        pair_commutator=pair_commutator(specification, include_empty=False),
        activity_reducer=graph_model.activity_polynomial,
        event_partition_certificate=graph_model.event_partition_certificate(),
        input_definition_sha256=contract.digest,
        dangling_certificate=graph_model.dangling_certificate(),
    )


def polynomial_jump_model(
    selected_queries: tuple[str, ...] = SELECTED_QUERIES,
) -> PolynomialJumpModel:
    """Return the authoritative composition-derived polynomial jump model."""

    return composition_compilation(selected_queries).model


def linear_closure(selected_queries=SELECTED_QUERIES, *, linear_limit=24):
    """Algorithm 1 from the molecularly derived typed rules and count reducer."""
    from mfobs.api import RatedRule, Observable
    from mfobs.core.counting_linear_closure import discover_linear_observables
    contract = analysis_contract(selected_queries)
    graph_model = BrusselatorGrowthGraphModel(input_definition_sha256=contract.digest)
    rules = tuple(RatedRule(r.rule, r.channel, identifier=r.identifier)
                  for r in graph_model.reaction_classes())
    expressions = graph_model.observable_expressions()
    queries = {name: Observable(name, terms=tuple(
        (c, expressions[name].representatives[k]) for k,c in expressions[name].coefficients.items()))
        for name in selected_queries}
    return discover_linear_observables(rules, queries, counter=graph_model,
        maximum_dimension=linear_limit, input_definition_sha256=contract.digest)


def analyze(
    *, linear_limit: int = 24,
    selected_queries: tuple[str, ...] = SELECTED_QUERIES,
    state_space=None,
    run_linear=True,
) -> dict:
    """Run both exact branches from the same normalized input contract."""

    from mfobs.models.brusselator_growth_chemistry import state_family as brusselator
    if state_space is not None:
        state_space.require(brusselator())
    contract = analysis_contract(selected_queries)
    linear = linear_closure(selected_queries, linear_limit=linear_limit) if run_linear else None
    compilation = composition_compilation(selected_queries)
    model = compilation.model
    pjc = analyze_polynomial_jump_closure(model)
    assert_same_input(contract, (linear, pjc) if linear is not None else (pjc,))
    projection = BrusselatorGrowthGraphModel(
        input_definition_sha256=contract.digest
    ).projection_record()
    return {
        "contract": contract.record(),
        "composition": {
            **compilation.evidence,
            "projection": projection,
        },
        "linear": linear,
        "pjc": pjc,
    }
