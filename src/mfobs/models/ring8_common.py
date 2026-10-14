"""Shared-contract linear and pJC analysis for the fixed Ring-8 model.

The linear branch is the existing conditional-composition worklist, seeded by
the single occupied-slot query ``P0``.  Its generated identity rules form an
exact coordinate basis on the 30 dihedral classes.  Independently, the pJC
branch applies the three literal conditional-DPO rules to one graph
representative of every class and lets the generic finite-language backend
derive all coordinate jumps, event propensities and the EMGF operator.

No transition matrix, jump vector, event partition or published equation is
supplied to this adapter.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from typing import Mapping

from mfobs.core.analysis_contract import (
    AnalysisContract,
    NormalizedRule,
    ObservableDeclaration,
    StateLanguageContract,
    assert_same_input,
)
from mfobs.core.finite_language_pjc import compile_finite_language_jump_model
from mfobs.core.finite_rule_projection import (
    FiniteRuleChannel,
    project_finite_rule_generators,
)
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
from mfobs.core.gml import read_rule_record
from mfobs.core.graph import has_edge
from mfobs.inputs import ModelManifest, load
from mfobs.models import ring_language as ring


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MODEL_PATH = PROJECT_ROOT / "examples" / "ring8" / "model.json"


def analysis_contract(inputs: ModelManifest | None = None) -> AnalysisContract:
    inputs = load(MODEL_PATH) if inputs is None else inputs
    rules = tuple(
        NormalizedRule(
            identifier=record["id"],
            rate=record["rate"],
            semantics="conditional-DPO",
            gml=inputs.file(record["file"]).read_text(encoding="utf-8"),
            conditions=(
                "NoEdge forbids a pre-existing bound edge at the target slot",
            ) if record["id"] != "unbind" else (),
        )
        for record in inputs.data["rules"]
    )
    selected = tuple(inputs.data["initial_observables"])
    observable_records = {
        record["id"]: record for record in inputs.data["observables"]
    }
    language_data = inputs.read("state_language")
    reduction_data = inputs.read("reduction")
    language = StateLanguageContract(
        identifier="ring-eight-fixed-sites-v1",
        carrier=(
            "one undirected eight-site C/P/L scaffold; only the eight optional "
            "bound edges vary"
        ),
        infinite=False,
        invariants=(
            "exactly eight C/P/L site triples and one next-cycle",
            "each site has exactly one cp and one cl edge",
            "each P/L pair has either zero or one bound edge",
            "rules and observables are invariant under the dihedral action",
        ),
        overlap_constraints=(
            "injective typed matches only",
            "transport and enforce every rule and observable NoEdge condition",
            "preserve the fixed scaffold and simple-graph constraints",
            "retain distinct cooperative-neighbour matches as multiplicities",
        ),
        preservation_argument=(
            "Each rule only inserts or removes one bound edge on the fixed "
            "scaffold, so all products remain in the declared language."
        ),
        semantic_quotient=(
            "the complete 30-class dihedral quotient of all 256 labelled bit states"
        ),
        finite_enumeration_role=(
            "complete proof domain for this finite model, not a sampled validation set"
        ),
    )
    return AnalysisContract(
        model_id="ring8",
        title="Ring-8 exhaustive correctness benchmark",
        rules=rules,
        state_language=language,
        observables=tuple(
            ObservableDeclaration(
                name,
                (
                    "number of occupied slots"
                    if name == "P0" else
                    "raw injective identity-pattern count on the fixed ring"
                ),
                "selected",
                (
                    "the sole identity-rule query in the checked-in model"
                    if name == "P0" else
                    "supplied identity-rule query"
                ),
                inputs.file(observable_records[name]["file"]).read_text(encoding="utf-8"),
            )
            for name in selected
        ),
        linear_seeds=selected,
        pjc_targets=selected,
        normalization={
            **{record["id"]: f'{record["per_match"]} per injective match'
               for record in inputs.data["rules"]},
            "coordinate_lift": reduction_data["basis"],
            "finite_language": language_data["assumption"],
        },
    )


def composition_linear_result(
    *, precomputed: Mapping, inputs: ModelManifest | None = None
) -> dict:
    """Validate the supplied composition result for the selected input."""

    inputs = load(MODEL_PATH) if inputs is None else inputs
    result = dict(precomputed)
    if result.get("input_definition_sha256") != inputs.definition_digest():
        raise ValueError("Ring-8 composition result does not match its manifest")
    selected = list(inputs.data["initial_observables"])
    if result.get("initial_observables") != selected:
        raise ValueError("Ring-8 discovery did not start from the selected queries")
    if not result.get("closed") or result.get("unknown") or result.get("pending"):
        raise ValueError("Ring-8 composition worklist is not exactly closed")
    dimension = len(result.get("basis", ()))
    if not dimension or result.get("basis_evaluation_rank") != dimension:
        raise ValueError("Ring-8 result is not an independent quotient basis")
    if dimension > 30:
        raise ValueError("Ring-8 basis exceeds the complete 30-state quotient dimension")
    return result


def quotient_states(result: Mapping) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(bit) for bit in state) for state in result["orbits"])


def finite_rule_projection(result: Mapping, *, inputs: ModelManifest | None = None):
    """Derive every rate-separated quotient transition by literal application."""

    inputs = load(MODEL_PATH) if inputs is None else inputs
    states = quotient_states(result)
    channels = tuple(
        FiniteRuleChannel(
            record["id"],
            read_rule_record(inputs.file(record["file"]).read_text(encoding="utf-8")),
            Fraction(record["per_match"]),
        )
        for record in inputs.data["rules"]
    )

    def representative(state: tuple[int, ...], _channel: str):
        return ring.slot_ring_graph(state)

    def project_target(_source: tuple[int, ...], _channel: str, target):
        bits = tuple(
            int(has_edge(target, 3 * site + 1, 3 * site + 2, "bound"))
            for site in range(8)
        )
        return ring.canonical_orbit(bits)

    return project_finite_rule_generators(
        states=states,
        channels=channels,
        representative=representative,
        project_target=project_target,
    )


def _observable_vectors(result: Mapping) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(value) for value in result["basis_vectors"][name])
        for name in result["basis"]
    )


def _validate_first_moments(result: Mapping, vectors, generators) -> dict:
    names = tuple(result["basis"])
    state_count = len(vectors[0])
    checks = 0
    for observable_index, observable in enumerate(names):
        for channel, generator in generators.items():
            expected = tuple(
                sum((generator[state][target] * vectors[observable_index][target]
                     for target in range(state_count)), Fraction())
                for state in range(state_count)
            )
            coefficients = result["equations_by_reaction"][observable][channel]
            actual = tuple(
                sum((Fraction(coefficients.get(target_name, "0"))
                     * vectors[target_index][state]
                     for target_index, target_name in enumerate(names)), Fraction())
                for state in range(state_count)
            )
            if actual != expected:
                raise AssertionError(f"Ring-8 first-moment mismatch: {observable}/{channel}")
            checks += state_count
    return {
        "status": "PASS",
        "exact_scalar_checks": checks,
        "comparison": (
            "literal conditional-DPO finite projection equals the independently "
            "composition-derived linear generator on all 30 quotient states"
        ),
    }


def analyze(
    *, composition_result: Mapping, inputs: ModelManifest | None = None
) -> dict:
    """Run both branches from one normalized input and cross-check their generators."""

    inputs = load(MODEL_PATH) if inputs is None else inputs
    contract = analysis_contract(inputs)
    result = composition_linear_result(precomputed=composition_result, inputs=inputs)
    states = quotient_states(result)
    vectors = _observable_vectors(result)
    projection = finite_rule_projection(result, inputs=inputs)
    compilation = compile_finite_language_jump_model(
        name=contract.model_id,
        states=states,
        observable_names=tuple(result["basis"]),
        observable_vectors=vectors,
        generators=projection.generators,
        input_definition_sha256=contract.digest,
        quotient_certificate={
            "method": "complete 30-class dihedral quotient",
            "all_256_labelled_states_covered": True,
            "finite_enumeration_role": contract.state_language.finite_enumeration_role,
        },
        source_rules={channel: channel for channel in projection.generators},
    )
    pjc = analyze_polynomial_jump_closure(compilation.model)
    pjc.update({
        "scope": "exact finite-language theorem on the complete Ring-8 quotient",
        "selected_targets": list(contract.pjc_targets),
        "coordinate_lift": {
            "selected": list(contract.pjc_targets),
            "generated_support": [name for name in result["basis"]
                                  if name not in contract.pjc_targets],
            "origin": (
                "conditional-composition linear discovery from P0; no support "
                "coordinate is supplied in the model manifest"
            ),
        },
        "finite_language_compilation": compilation.evidence,
    })
    linear = {
        "schema_version": 1,
        "status": "CLOSED",
        "analysis": "linear-observable-discovery",
        "input_definition_sha256": contract.digest,
        "seeds": list(contract.linear_seeds),
        "basis": list(result["basis"]),
        "dimension": len(result["basis"]),
        "scope": "exact on the complete finite Ring-8 language",
        "composition_artifact": "explicit composition result",
        "ordered_composition_calls": result["ordered_composition_calls"],
        "raw_compositions": result["raw_compositions"],
    }
    assert_same_input(contract, (linear, pjc))
    return {
        "schema_version": 1,
        "contract": contract.record(),
        "composition": {
            "schema_version": 1,
            "status": "PASS",
            "jump_projection": "DPO",
            "scope": "linear first-commutator observable discovery",
            "full_event_partition_composition_derived": False,
            "generated_support_coordinates": pjc["coordinate_lift"]["generated_support"],
        },
        "finite_rule_projection": projection.evidence,
        "cross_branch_validation": _validate_first_moments(
            result, vectors, projection.generators
        ),
        "linear": linear,
        "pjc": pjc,
    }
