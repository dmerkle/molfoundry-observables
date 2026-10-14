"""Compile finite polynomial jump models from graph-rule compositions.

Applications provide a finite, certified partition of rule matches, diagonal
graph observables, and an exact reducer for each partition's activity identity.
This module derives every observable jump by rule composition.  It accepts a
reaction class only when the jump-closed commutator is an exact integer
multiple of that class's activity identity under the selected DPO or SqPO
jump projection.  Thus event deltas are outputs, not application-supplied
stoichiometric tables.

The compiler is independent of chemistry, words, trees, and observable names.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Callable, Mapping

from mfobs.core.engine import dpo_jump_closure
from mfobs.core.jump_closure_worklist import (
    RuleExpression,
    bilinear_commutator,
    jump_close_expression,
)
from mfobs.core.identity import certified_input_identity, identity_contradictions
from mfobs.core.observables import parse_observable
from mfobs.core.polynomial_jump_closure import (
    JumpEvent,
    JumpModelProvenance,
    Polynomial,
    PolynomialJumpModel,
    PolynomialRing,
)


class NonconstantCompositionJump(ValueError):
    """A declared match class does not have one constant observable jump."""


@dataclass(frozen=True)
class GraphReactionClass:
    """One graph-defined, constant-jump partition of a normalized rule."""

    identifier: str
    channel: str
    rule: object
    source_rule: str
    predicate: str
    metadata: Mapping


@dataclass(frozen=True)
class CompositionJumpCompilation:
    """The compiled jump model and its inspectable composition evidence."""

    model: PolynomialJumpModel
    evidence: Mapping


def _digest(key: str) -> str:
    return sha256(key.encode()).hexdigest()


def _expression_record(expression: RuleExpression) -> dict:
    """Serialize coefficients together with the retained raw rule witnesses."""

    record = expression.record()
    witnesses = {
        _digest(key): expression.representatives[key].getGMLString()
        for key in expression.coefficients
    }
    return {
        **record,
        "terms": [
            {**term, "rule_gml": witnesses[term["term_id"]]}
            for term in record["terms"]
        ],
    }


def compile_composition_jump_model(
    *,
    name: str,
    ring: PolynomialRing,
    observables: Mapping[str, RuleExpression],
    reactions: tuple[GraphReactionClass, ...],
    pair_commutator: Callable[[object, object], object],
    activity_reducer: Callable[[PolynomialRing, GraphReactionClass, object], Polynomial],
    event_partition_certificate: Mapping,
    input_definition_sha256: str = "",
    jump_closure: Callable[[object], object] = dpo_jump_closure,
    jump_projection_name: str = "DPO",
    composition_backend_name: str = "MolFoundry DPO common-overlap composition",
    dangling_certificate: Callable[[object, tuple[int, ...]], None] | None = None,
) -> CompositionJumpCompilation:
    """Derive all constant jumps and construct a polynomial jump model.

    ``activity_reducer`` is the sole state-language-specific operation: it
    maps the selected activity identity of a graph rule to its exact
    polynomial match count.  Observable changes are obtained only from
    compositions.  ``jump_closure`` selects DPO or SqPO semantics explicitly;
    it is used consistently for activities and observable commutators.
    A DPO jump closure that deletes a vertex counts matches only under its
    dangling condition; it is reduced to an input count only when
    ``dangling_certificate(graph, deleted)`` proves that condition on the
    declared language (see :func:`certified_input_identity`).
    """

    if tuple(observables) != ring.variables:
        raise ValueError("observable order must equal the polynomial ring variables")
    for name_of_observable, expression in observables.items():
        if not expression.coefficients or any(
            not parse_observable(rule).jump_identity
            for rule in expression.representatives.values()
        ):
            raise ValueError(
                f"{name_of_observable}: coordinates must be nonzero identity-rule expressions"
            )
    if not reactions or len({reaction.identifier for reaction in reactions}) != len(reactions):
        raise ValueError("reaction classes must be nonempty and uniquely named")

    events = []
    records = []
    for reaction in reactions:
        reaction_expression = RuleExpression.from_rules((reaction.rule,))
        activity = jump_closure(
            reaction.rule,
            name=f"activity_{reaction.identifier}",
        )
        normalized_activity = certified_input_identity(activity, dangling_certificate)
        activity_key = normalized_activity.canonKey
        activity_zero_reasons = identity_contradictions(normalized_activity)
        if activity_zero_reasons:
            raise ValueError(
                f"{reaction.identifier}: contradictory activity identity"
            )
        propensity = activity_reducer(ring, reaction, normalized_activity)
        if propensity.is_zero:
            raise ValueError(f"{reaction.identifier}: activity reduces to zero")

        delta = []
        observable_records = {}
        for observable_name, observable_expression in observables.items():
            raw, pair_records = bilinear_commutator(
                observable_expression,
                reaction_expression,
                pair_commutator,
            )
            jumped, jump_trace = jump_close_expression(
                raw,
                jump_closure=jump_closure,
            )
            eliminated = {}
            for key in jumped.coefficients:
                reasons = identity_contradictions(jumped.representatives[key])
                if reasons:
                    eliminated[key] = reasons
            if eliminated:
                jumped = RuleExpression(
                    {
                        key: coefficient
                        for key, coefficient in jumped.coefficients.items()
                        if key not in eliminated
                    },
                    {
                        key: rule
                        for key, rule in jumped.representatives.items()
                        if key not in eliminated
                    },
                )
            normalized_coefficients = {}
            normalized_representatives = {}
            for key, coefficient in jumped.coefficients.items():
                normalized = certified_input_identity(
                    jumped.representatives[key], dangling_certificate
                )
                normalized_key = normalized.canonKey
                normalized_coefficients[normalized_key] = (
                    normalized_coefficients.get(normalized_key, 0) + coefficient
                )
                normalized_representatives.setdefault(normalized_key, normalized)
            jumped = RuleExpression(
                normalized_coefficients,
                normalized_representatives,
            )
            unexpected = {
                key: coefficient
                for key, coefficient in jumped.coefficients.items()
                if key != activity_key and coefficient
            }
            if unexpected:
                raise NonconstantCompositionJump(
                    f"{reaction.identifier}/{observable_name}: jump-closed "
                    f"commutator has {len(unexpected)} non-activity terms"
                )
            coefficient = jumped.coefficients.get(activity_key, 0)
            if getattr(coefficient, "denominator", 1) != 1:
                raise NonconstantCompositionJump(
                    f"{reaction.identifier}/{observable_name}: non-integral jump"
                )
            value = int(coefficient)
            delta.append(value)
            observable_records[observable_name] = {
                "delta": value,
                "raw_commutator": _expression_record(raw),
                "jump_closed": _expression_record(jumped),
                "composition_pairs": pair_records,
                "jump_closure_trace": jump_trace,
                "contradictory_terms_eliminated": {
                    _digest(key): list(reasons)
                    for key, reasons in eliminated.items()
                },
            }

        events.append(JumpEvent(
            reaction.identifier,
            reaction.channel,
            tuple(delta),
            propensity,
            reaction.source_rule,
            reaction.predicate,
        ))
        records.append({
            "event_class": reaction.identifier,
            "source_rule": reaction.source_rule,
            "rate_channel": reaction.channel,
            "predicate": reaction.predicate,
            "partition_metadata": dict(reaction.metadata),
            "structural_rule_gml": reaction.rule.getGMLString(),
            "activity_identity_sha256": _digest(activity_key),
            "activity_identity_gml": normalized_activity.getGMLString(),
            "activity_polynomial": ring.record(propensity),
            "composition_derived_delta": {
                variable: value for variable, value in zip(ring.variables, delta)
                if value
            },
            "observable_commutators": observable_records,
        })

    certificate = {
        **dict(event_partition_certificate),
        "jump_derivation": (
            "Every coordinate jump is the exact coefficient of the reaction "
            "class activity identity in its selected jump-closed commutator."
        ),
        "composition_used": True,
        "supplied_delta_table": False,
    }
    model = PolynomialJumpModel(
        name,
        ring,
        tuple(events),
        certificate,
        input_definition_sha256,
        JumpModelProvenance(
            backend=composition_backend_name,
            event_partition_origin=(
                "application-supplied all-state structural match partition"
            ),
            jump_origin=(
                "exact activity-identity coefficient of each composed "
                "observable commutator"
            ),
            activity_origin="application-supplied exact state-language reducer",
            jump_projection=jump_projection_name,
            composition_scope="all event-coordinate jumps",
            all_jumps_composition_derived=True,
            supplied_delta_table=False,
        ),
    )
    return CompositionJumpCompilation(
        model,
        {
            "schema_version": 2,
            "status": "PASS",
            "backend": composition_backend_name,
            "jump_projection": jump_projection_name,
            "arithmetic": "exact integral activity-identity coefficients",
            "input_definition_sha256": input_definition_sha256,
            "observable_coordinates": list(ring.variables),
            "event_classes": records,
        },
    )


def compare_polynomial_jump_models(
    derived: PolynomialJumpModel,
    oracle: PolynomialJumpModel,
) -> dict:
    """Compare two jump models without using either one to construct the other."""

    failures = []
    def ring_signature(value: PolynomialRing):
        return (
            value.variables,
            value.boolean_indices,
            value.exclusive_groups,
            tuple(sorted(value.boolean_products.items())),
        )

    if ring_signature(derived.ring) != ring_signature(oracle.ring):
        failures.append({
            "kind": "coordinate-algebra",
            "derived": list(derived.ring.variables),
            "oracle": list(oracle.ring.variables),
        })
    derived_events = {event.identifier: event for event in derived.events}
    oracle_events = {event.identifier: event for event in oracle.events}
    if set(derived_events) != set(oracle_events):
        failures.append({
            "kind": "event-identifiers",
            "derived": sorted(derived_events),
            "oracle": sorted(oracle_events),
        })
    for identifier in sorted(set(derived_events) & set(oracle_events)):
        left, right = derived_events[identifier], oracle_events[identifier]
        comparisons = {
            "channel": (left.channel, right.channel),
            "delta": (left.delta, right.delta),
            "propensity": (left.propensity.terms, right.propensity.terms),
            "source_rule": (left.rule, right.rule),
        }
        for field, (actual, expected) in comparisons.items():
            if actual != expected:
                failures.append({
                    "kind": field,
                    "event": identifier,
                    "derived": str(actual),
                    "oracle": str(expected),
                })
    return {
        "status": "PASS" if not failures else "FAIL",
        "scope": (
            "Independent equality check: composition-derived graph model "
            "versus semantic oracle"
        ),
        "event_classes": len(derived_events),
        "failures": failures,
    }
