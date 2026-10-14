"""Algorithm 2: extract finite events from closed polynomial support.

The backend justifies uniform source-rule matches and supplies a finite,
nonempty, disjoint, exhaustive partition of enabled coordinate states. Each
class has a polynomial indicator and substitutions making every jump constant.
These are obligations on the declared graph family, not consequences of a
finite test. This module knows no graph labels, constructors or pattern kinds.
"""
from dataclasses import dataclass
from typing import Mapping

import sympy as sp

from mfobs.core.pattern_support import UnsupportedReduction, polynomial_record


@dataclass(frozen=True)
class EventClass:
    """Backend-provided class; substitutions are identities within that class."""
    substitutions: Mapping
    indicator: sp.Expr
    predicate: Mapping


def constant_partition(source, jumps):
    """One enabled class, supported only when all derived jumps are constant."""
    if any(value.free_symbols for value in jumps.values()):
        raise UnsupportedReduction("backend cannot partition state-dependent jumps")
    return (EventClass({}, sp.Integer(1), {}),)


def extract_uniform_events(support, sources, *, partition, scope,
                           uniform_match_argument, partition_argument):
    """Compute F/a, refine by backend classes, and emit integer jumps and rates.

    A backend must justify its arguments on all enabled states. They are
    recorded declarations, not machine-checked proofs. Failure returns UNKNOWN
    with no partial event model advertised as complete.
    """
    if support.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status": "UNKNOWN_SUPPORT", "events": []}
    if not uniform_match_argument or not partition_argument:
        return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "events": []}
    events, obligations = [], []
    try:
        for source in sources:
            activity = sp.sympify(support.activities[source.identifier])
            if activity == 0:
                continue
            jumps = {name: sp.cancel(support.images[name][source.identifier] / activity)
                     for name in support.patterns}
            for value in jumps.values():
                polynomial_record(value)
            classes = list(partition(source, jumps))
            if not classes:
                raise UnsupportedReduction("empty partition for a nonzero activity")
            for cell in classes:
                polynomial_record(cell.indicator)
                delta = {name: sp.cancel(value.subs(cell.substitutions, simultaneous=True))
                         for name, value in jumps.items()}
                if any(value.is_Integer is not True for value in delta.values()):
                    raise UnsupportedReduction("partition does not give constant integral jumps")
                if not any(delta.values()):
                    continue
                events.append({"source_rule": source.identifier, "rate": source.rate,
                               "predicate": dict(cell.predicate),
                               "jump": {name: int(value) for name, value in delta.items() if value},
                               "activity": polynomial_record(activity * cell.indicator)})
            obligations.append({"source_rule": source.identifier,
                "uniform_match_argument": uniform_match_argument,
                "jump_polynomials": {name: str(value) for name, value in jumps.items()}})
        return {"status": "FINITE_EVENTS", "events": events,
                "uniform_match_obligations": obligations, "scope": scope,
                "partition_argument": partition_argument,
                "event_partition_origin": partition_argument}
    except (UnsupportedReduction, sp.PolynomialError) as exc:
        return {"status": "UNKNOWN_EVENT_EXTRACTION", "reason": str(exc), "events": []}


@dataclass(frozen=True)
class UniformEventCompilation:
    model: object
    evidence: Mapping


def compile_event_model(extracted, *, name, ring, certificate,
                        input_definition_sha256="", provenance):
    """Serialize extracted events in the common numerical/EMGF representation."""
    from mfobs.core.polynomial_jump_closure import JumpEvent, PolynomialJumpModel
    if extracted["status"] != "FINITE_EVENTS":
        raise UnsupportedReduction(extracted.get("reason", extracted["status"]))
    events = []
    for index, event in enumerate(extracted["events"]):
        propensity = ring.add(*[
            ring.scale(ring.monomial(term["powers"]), term["coefficient"])
            for term in event["activity"]["terms"]])
        if propensity.is_zero:
            continue
        events.append(JumpEvent(
            f"{event['source_rule']}-support-{index}", event["rate"],
            tuple(event["jump"].get(n, 0) for n in ring.variables), propensity,
            event["source_rule"],
            " and ".join(f"{n}={v}" for n, v in event["predicate"].items()) or "true"))
    model = PolynomialJumpModel(name, ring, tuple(events), certificate,
                                input_definition_sha256, provenance)
    return UniformEventCompilation(model, extracted)
