"""Shared, inspectable contract for observable analyses.

Every model-facing analysis starts from the same four ingredients: literal
rules, a state language, structural obligations on that language, and named
seed observables.  The two consumers are deliberately separate:

* linear discovery studies first generator images and may return ``UNKNOWN``;
* polynomial jump closure constructs an EMGF operator from constant jumps and
  polynomial event propensities.

This module contains metadata and validation only.  It neither knows a graph
model nor chooses observables on behalf of an example.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from typing import Mapping, Sequence


def _nonempty(value: str, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")
    return value


@dataclass(frozen=True)
class StateLanguageContract:
    """All-size semantic domain used by an exact symbolic reducer."""

    identifier: str
    carrier: str
    infinite: bool
    invariants: tuple[str, ...]
    overlap_constraints: tuple[str, ...]
    preservation_argument: str
    semantic_quotient: str
    finite_enumeration_role: str = "validation only"

    def __post_init__(self) -> None:
        for field in (
            "identifier", "carrier", "preservation_argument",
            "semantic_quotient", "finite_enumeration_role",
        ):
            _nonempty(getattr(self, field), field)
        if not self.invariants:
            raise ValueError("a state language needs at least one invariant")
        if not self.overlap_constraints:
            raise ValueError("a state language needs explicit overlap constraints")

    def record(self) -> dict:
        return {
            "id": self.identifier,
            "carrier": self.carrier,
            "infinite": self.infinite,
            "invariants": list(self.invariants),
            "overlap_constraints": list(self.overlap_constraints),
            "preservation_argument": self.preservation_argument,
            "semantic_quotient": self.semantic_quotient,
            "finite_enumeration_role": self.finite_enumeration_role,
        }


@dataclass(frozen=True)
class NormalizedRule:
    """One rate-separated rule in the common analysis input."""

    identifier: str
    rate: str
    semantics: str
    gml: str
    conditions: tuple[str, ...] = ()
    boundary: bool = False

    def __post_init__(self) -> None:
        _nonempty(self.identifier, "rule identifier")
        _nonempty(self.rate, "rate")
        if self.semantics != "conditional-DPO":
            raise ValueError("normalized examples currently require conditional-DPO semantics")
        if "rule [" not in self.gml or "ruleID" not in self.gml:
            raise ValueError(f"{self.identifier}: literal GML rule is missing")

    def record(self) -> dict:
        return {
            "id": self.identifier,
            "rate": self.rate,
            "semantics": self.semantics,
            "conditions": list(self.conditions),
            "boundary_operation": self.boundary,
            "gml_sha256": sha256(self.gml.encode()).hexdigest(),
        }


@dataclass(frozen=True)
class ObservableDeclaration:
    """A query, support coordinate, or generated diagonal observable."""

    name: str
    meaning: str
    role: str
    origin: str
    gml: str | None = None

    def __post_init__(self) -> None:
        _nonempty(self.name, "observable name")
        _nonempty(self.meaning, "observable meaning")
        if self.role not in {"selected", "support", "auxiliary", "derived"}:
            raise ValueError(f"{self.name}: unsupported observable role {self.role!r}")
        if self.gml is not None and "rule [" not in self.gml:
            raise ValueError(f"{self.name}: observable GML is not an identity rule")

    def record(self) -> dict:
        result = {
            "name": self.name,
            "meaning": self.meaning,
            "role": self.role,
            "origin": self.origin,
        }
        if self.gml is not None:
            result["gml_sha256"] = sha256(self.gml.encode()).hexdigest()
        return result


@dataclass(frozen=True)
class AnalysisContract:
    """Single normalized input consumed by both analysis branches."""

    model_id: str
    title: str
    rules: tuple[NormalizedRule, ...]
    state_language: StateLanguageContract
    observables: tuple[ObservableDeclaration, ...]
    linear_seeds: tuple[str, ...]
    pjc_targets: tuple[str, ...]
    normalization: Mapping[str, str]

    def __post_init__(self) -> None:
        _nonempty(self.model_id, "model_id")
        _nonempty(self.title, "title")
        rule_ids = [rule.identifier for rule in self.rules]
        observable_names = [observable.name for observable in self.observables]
        if not rule_ids or len(rule_ids) != len(set(rule_ids)):
            raise ValueError("rule identifiers must be nonempty and unique")
        if not observable_names or len(observable_names) != len(set(observable_names)):
            raise ValueError("observable names must be nonempty and unique")
        known = set(observable_names)
        if not self.linear_seeds or set(self.linear_seeds) - known:
            raise ValueError("linear seeds must name declared observables")
        if not self.pjc_targets or set(self.pjc_targets) - known:
            raise ValueError("pJC targets must name declared observables")

    @property
    def digest(self) -> str:
        encoded = json.dumps(self.record(include_digest=False), sort_keys=True).encode()
        return sha256(encoded).hexdigest()

    def record(self, *, include_digest: bool = True) -> dict:
        result = {
            "schema_version": 1,
            "model": self.model_id,
            "title": self.title,
            "normalized_rules": [rule.record() for rule in self.rules],
            "state_language": self.state_language.record(),
            "observables": [observable.record() for observable in self.observables],
            "analysis_branches": {
                "linear_observable_discovery": {
                    "seeds": list(self.linear_seeds),
                    "operation": "first commutator followed by diagonal jump projection",
                },
                "polynomial_jump_closure_emgf": {
                    "targets": list(self.pjc_targets),
                    "operation": "all nested diagonal commutators via polynomial event propensities",
                },
            },
            "normalization": dict(self.normalization),
        }
        if include_digest:
            result["definition_sha256"] = self.digest
        return result


def assert_same_input(contract: AnalysisContract, analyses: Sequence[Mapping]) -> None:
    """Fail if branch artifacts were not derived from the same contract."""

    for result in analyses:
        if result.get("input_definition_sha256") != contract.digest:
            raise ValueError("analysis result does not match the normalized input contract")
