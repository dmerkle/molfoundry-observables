"""Native rule/count declarations and the public analysis dispatcher."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Mapping
from typing import Protocol

import molfoundry as mf
from mfobs.state_space import StateSpace


@dataclass(frozen=True)
class RatedRule:
    """A MolFoundry rule with a symbolic rate channel."""

    rule: mf.Rule
    rate: str = "k"
    per_match: Fraction = Fraction(1)
    identifier: str | None = None

    def __post_init__(self):
        factor = Fraction(self.per_match)
        if factor <= 0:
            raise ValueError("per_match must be positive")
        object.__setattr__(self, "per_match", factor)
        if self.identifier is None:
            object.__setattr__(self, "identifier", self.rule.name)

    @property
    def name(self) -> str:
        return self.rule.name


@dataclass(frozen=True)
class Observable:
    """A named graph count, conditional identity count, or sum of such counts."""

    name: str
    pattern: mf.Graph | mf.Rule | None = None
    terms: tuple[tuple[Fraction, mf.Rule], ...] = ()

    def identity_terms(self) -> tuple[tuple[Fraction, mf.Rule], ...]:
        from mfobs.core.rules import is_identity

        if self.pattern is not None:
            if self.terms:
                raise ValueError("supply a pattern or identity terms, not both")
            rule = (self.pattern if isinstance(self.pattern, mf.Rule) else
                    mf.Rule.identity(self.pattern, self.name, add=False))
            terms = ((Fraction(1), rule),)
        else:
            terms = tuple((Fraction(coefficient), rule) for coefficient, rule in self.terms)
        if not terms or any(not c or not is_identity(r) for c, r in terms):
            raise ValueError(f"{self.name}: nonzero identity-rule terms required")
        return terms


@dataclass(frozen=True)
class ModelInputs:
    """The supplied Ω, rated rules and selected observable expressions.

    Packaged families have specialised counting backends. Changes outside a
    backend's supported input are rejected before analysis, never ignored.
    """

    omega: StateSpace
    rules: Mapping[str, RatedRule]
    observables: Mapping[str, Observable]

    def __post_init__(self):
        if not isinstance(self.omega, StateSpace):
            raise TypeError("omega must be a StateSpace declaration")
        if not self.rules or not self.observables:
            raise ValueError("rules and selected observables must be nonempty")
        object.__setattr__(self, "rules", MappingProxyType(dict(self.rules)))
        object.__setattr__(self, "observables", MappingProxyType(dict(self.observables)))
        for name, query in self.observables.items():
            if name != query.name:
                raise ValueError("observable key must equal its name")
            query.identity_terms()

    def rule_gmls(self) -> dict[str, str]:
        return {name: value.rule.getGMLString() for name, value in self.rules.items()}

    def query_gmls(self) -> dict[str, str]:
        result = {}
        for name, query in self.observables.items():
            terms = query.identity_terms()
            for index, (_, rule) in enumerate(terms):
                key = name if len(terms) == 1 else f"{name}:{index + 1}"
                result[key] = rule.getGMLString()
        return result

    def signature(self):
        def expression(query):
            total = {}
            for coefficient, rule in query.identity_terms():
                total[rule.canonKey] = total.get(rule.canonKey, Fraction(0)) + coefficient
            return {key: value for key, value in total.items() if value}
        return (
            self.omega.record(),
            {name: (r.rule.canonKey, r.rate, r.per_match, r.identifier) for name, r in self.rules.items()},
            {name: expression(query) for name, query in self.observables.items()},
        )

    def require(self, expected: "ModelInputs") -> None:
        """Check the input scope justified by a specialised counting backend."""
        if self.signature() != expected.signature():
            raise ValueError("unsupported Ω, rules or observables for this counting backend")

    def validate(self) -> None:
        for declaration in self.rules.values():
            self.omega.validate_rule(declaration.rule)
        for query in self.observables.values():
            for _, identity in query.identity_terms():
                self.omega.validate_rule(identity)


@dataclass(frozen=True)
class JumpClosureResult:
    """Algorithm 2: discovered support, event model and pJC status."""

    support: object
    compilation: object
    pjc: dict


@dataclass(frozen=True)
class AnalysisBranch:
    """One computed outcome, including its status and inspectable evidence."""

    data: Mapping

    @property
    def status(self) -> str:
        return self.data["status"]

    def record(self) -> dict:
        return dict(self.data)


@dataclass(frozen=True)
class AnalysisResult:
    """Common public result; an absent branch means it was not constructed.

    The branch fields are independent. Reporting a closed polynomial support
    does not imply linear closure or an event model. ``record()`` exports the
    model's numerical/evidence record for the existing rendering consumers.
    """

    inputs: ModelInputs
    linear: AnalysisBranch | None
    support: AnalysisBranch | None
    pjc: AnalysisBranch | None
    evidence: Mapping

    @classmethod
    def from_record(cls, inputs, record):
        def branch(value):
            return AnalysisBranch(MappingProxyType(dict(value))) if value is not None else None
        linear = record.get("linear")
        if linear is None and "basis" in record and "status" in record:
            linear = record
        support = record.get("support_discovery") or record.get("synthesis", {}).get("support")
        return cls(inputs, branch(linear), branch(support), branch(record.get("pjc")),
                   MappingProxyType(dict(record)))

    def record(self) -> dict:
        return dict(self.evidence)

    def show(self) -> None:
        for name, result in (("Linear closure", self.linear), ("Polynomial support", self.support),
                             ("Polynomial jump closure", self.pjc)):
            print(name + ":", result.status if result is not None else "NOT_RUN")


class AnalysisBackend(Protocol):
    """Execute the common inputs using explicitly selected family semantics."""

    def run(self, inputs: ModelInputs, **options) -> dict:
        ...


def analyze(inputs: ModelInputs, *, backend: AnalysisBackend, **options) -> AnalysisResult:
    """Run the selected backend on Ω, native rules and selected observables."""
    if not isinstance(inputs, ModelInputs):
        raise TypeError("analyze requires ModelInputs and an explicit backend")
    if {"omega", "rules", "query", "queries", "observables", "state_space", "language"} & options.keys():
        raise TypeError("select Ω, rules and observables in ModelInputs")
    inputs.validate()
    return AnalysisResult.from_record(inputs, backend.run(inputs, **options))


__all__ = ["AnalysisBackend", "AnalysisBranch", "AnalysisResult", "JumpClosureResult",
           "Observable", "ModelInputs", "RatedRule", "analyze"]
