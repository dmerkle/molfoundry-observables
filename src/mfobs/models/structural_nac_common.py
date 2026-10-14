"""Shared analysis contract for the unbounded growing-ring oscillator.

This example deliberately reaches the boundary of the common architecture.
The linear branch has an exact symbolic local-word reducer but its FIFO
worklist remains open.  The selected one-site coordinates are not a Markov
quotient, and every finite word prefix seen so far grows a larger frontier.
Consequently the pJC branch returns an inspectable ``UNKNOWN`` certificate
instead of inventing a finite event table or applying the finite-language
backend outside its preconditions.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from itertools import product
from pathlib import Path
from typing import Mapping

from mfobs.core.analysis_contract import (
    AnalysisContract,
    NormalizedRule,
    ObservableDeclaration,
    StateLanguageContract,
    assert_same_input,
)
from mfobs.inputs import load
from mfobs.models import structural_nac_oscillator as chemistry


PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_ROOT = PROJECT_ROOT / "examples" / "structural-nac-oscillator"


def analysis_contract(inputs=None) -> AnalysisContract:
    inputs = load(EXAMPLE_ROOT / "model.json") if inputs is None else inputs
    rules = tuple(
        NormalizedRule(
            identifier=record["id"],
            rate=record["rate"],
            semantics="conditional-DPO",
            gml=inputs.file(record["file"]).read_text(encoding="utf-8"),
            conditions=(
                f"NoEdge requires the inspected neighbour not to carry closure e{chemistry.RULE_TRANSITIONS[record['id']][3]}"
            ) if record["id"] in chemistry.RULE_TRANSITIONS
                 and chemistry.RULE_TRANSITIONS[record["id"]][3] is not None else (),
            boundary=record["id"] in chemistry.GROW_STATES,
        )
        for record in inputs.data["rules"]
    )
    selected_names = tuple(inputs.data["initial_observables"])
    observable_files = {item["id"]: item["file"] for item in inputs.data["observables"]}
    observables = tuple(
        ObservableDeclaration(
            name,
            (
                f"number of sites carrying covalent closure topology {name}"
                if len(name) == 2 else
                "number of oriented consecutive occurrences of "
                + " → ".join(f"Q{value}" for value in name[1:])
            ),
            "selected",
            (
                "literal one-site identity rule from the model manifest"
                if len(name) == 2 else
                "literal connected-word identity selected through the common query interface"
            ),
            inputs.file(observable_files[name]).read_text(encoding="utf-8"),
        )
        for name in selected_names
    )
    language_data = inputs.read("state_language")
    return AnalysisContract(
        model_id="structural-nac-oscillator",
        title="Growing-ring structural NoEdge oscillator",
        rules=rules,
        state_language=StateLanguageContract(
            identifier="unbounded-oriented-covalent-switch-macrocycles-v1",
            carrier=(
                "one oriented covalent macrocycle of every size N >= 7 over "
                "the local closure alphabet Q0,Q1,Q2, with chemostatted capped precursors"
            ),
            infinite=True,
            invariants=tuple(language_data["axioms"]),
            overlap_constraints=(
                "injective typed conditional-DPO matches only",
                "transport all three active neighbouring-site NoEdge conditions",
                "preserve one oriented macrocycle and exactly one closure edge per site",
                "retain distinct local matches and growth templates as multiplicities",
                "reduce only rule/observable windows that fit in the declared ring",
            ),
            preservation_argument=(
                "Conversion and swap rules preserve the scaffold and ring size. "
                "Each growth rule consumes a complete capped precursor, releases "
                "mapped water and inserts exactly one valid repeat, hence N -> N+1."
            ),
            semantic_quotient=(
                "oriented cyclic Q0/Q1/Q2 words modulo rotation; connected word "
                "counts are the exact symbolic coordinate family"
            ),
            finite_enumeration_role=(
                "small rings validate generated rows only; they cannot certify "
                "closure of the unbounded all-size language"
            ),
        ),
        observables=observables,
        linear_seeds=selected_names,
        pjc_targets=selected_names,
        normalization={
            **{record["id"]: f'{record["per_match"]} per injective match'
               for record in inputs.data["rules"]},
            "word_counts": "one per oriented cyclic start",
        },
    )


def _selected_coordinate_witness() -> dict:
    """Find two same-count rings with different literal-rule activities."""

    groups = defaultdict(list)
    for state in product(range(3), repeat=7):
        signature = tuple(state.count(value) for value in range(3))
        activities = tuple(
            sum(1 for _ in chemistry.transition_events(state, reaction))
            for reaction in chemistry.REACTION_ORDER
        )
        for previous_state, previous_activities in groups[signature]:
            if previous_activities != activities:
                differing = [
                    reaction for reaction, left, right in zip(
                        chemistry.REACTION_ORDER, previous_activities, activities
                    ) if left != right
                ]
                return {
                    "status": "WITNESS",
                    "ring_size": 7,
                    "selected_coordinate_signature": {
                        "Q0": signature[0], "Q1": signature[1], "Q2": signature[2]
                    },
                    "state_a": list(previous_state),
                    "state_b": list(state),
                    "activities_a": dict(zip(chemistry.REACTION_ORDER, previous_activities)),
                    "activities_b": dict(zip(chemistry.REACTION_ORDER, activities)),
                    "differing_channels": differing,
                    "conclusion": (
                        "Q0,Q1,Q2 alone are not a strongly lumpable event quotient"
                    ),
                }
        groups[signature].append((state, activities))
    raise AssertionError("expected a selected-coordinate non-lumpability witness")


def analyze(
    *,
    linear_result: Mapping,
    profile: Mapping,
    inputs=None,
) -> dict:
    """Return the two branches with an explicit, evidence-backed pJC boundary."""

    inputs = load(EXAMPLE_ROOT / "model.json") if inputs is None else inputs
    contract = analysis_contract(inputs)
    result = dict(linear_result)
    profile_result = dict(profile)
    if result.get("input_definition_sha256") != inputs.definition_digest():
        raise ValueError("growing-ring discovery result does not match its manifest")
    if result.get("initial_observables") != list(contract.linear_seeds):
        raise ValueError("growing-ring discovery seeds differ from the selected queries")
    if result.get("closed") or not result.get("pending"):
        raise ValueError("growing-ring oracle unexpectedly claims finite closure")

    linear = {
        "schema_version": 1,
        "status": result.get('algorithm1_status', "UNKNOWN_EXPANSION_LIMIT"),
        "analysis": "linear-observable-discovery",
        "input_definition_sha256": contract.digest,
        "seeds": list(contract.linear_seeds),
        "basis": list(result["basis"]),
        "dimension": len(result["basis"]),
        "pending": list(result["pending"]),
        "obligations": list(result.get("obligations", [])),
        "basis_expressions": result.get("word_coordinates", {}),
        "scope": "exact generated prefix on the certified all-size word language",
        "expanded_coordinates": len(result["steps"]),
        "maximum_generated_word_length": result["maximum_generated_word_length"],
        "profile": {
            key: profile_result[key] for key in (
                "budget", "status", "basis_count", "pending_count",
                "maximum_word_length", "method", "interpretation",
            )
        },
    }
    witness = (
        _selected_coordinate_witness()
        if contract.linear_seeds == ("Q0", "Q1", "Q2") else None
    )
    pjc = {
        "schema_version": 1,
        "status": "UNKNOWN_UNBOUNDED_COORDINATE_LIFT",
        "analysis": "polynomial-jump-closure-emgf",
        "input_definition_sha256": contract.digest,
        "scope": "all oriented macrocycles of every size N >= 7",
        "selected_targets": list(contract.pjc_targets),
        "event_class_count": None,
        "operator_span": None,
        "operator": None,
        "finite_language_backend_used": False,
        "attempt_provenance": {
            "schema_version": 1,
            "backend": "not constructed: finite event quotient precondition failed",
            "event_partition_origin": "no finite all-size event partition certified",
            "jump_origin": "not derived beyond exact local-word generator prefixes",
            "activity_origin": "literal-rule activities expose selected-coordinate non-lumpability",
            "jump_projection": "DPO local-word semantics",
            "composition_scope": "exact first-generator word rows only",
            "all_jumps_composition_derived": False,
            "supplied_delta_table": False,
        },
        "reason": (
            "Exact support coordinates are connected words, and the common linear "
            "worklist continues to generate longer words with a nonempty frontier. "
            "No finite all-size event quotient has therefore been certified for the "
            "selected query family."
        ),
        "selected_coordinate_non_lumpability": witness,
        "prohibited_shortcut": (
            "Finite N enumeration or the approximate R1-R3 Markov closures are not "
            "accepted as an exact pJC certificate."
        ),
    }
    assert_same_input(contract, (linear, pjc))
    return {
        "schema_version": 1,
        "contract": contract.record(),
        "linear": linear,
        "pjc": pjc,
        "composition": {
            "schema_version": 1,
            "status": "MODEL_SPECIFIC_EXACT_REDUCER",
            "jump_projection": "DPO local-word semantics",
            "scope": "exact first-generator rows for connected word observables",
            "full_event_partition_composition_derived": False,
            "reason": (
                "Literal graph rules and guards are validated, while the required "
                "all-size quotient reduction is the model-specific cyclic-word theorem."
            ),
        },
        "source_manifest_sha256": sha256((EXAMPLE_ROOT / "model.json").read_bytes()).hexdigest(),
    }
