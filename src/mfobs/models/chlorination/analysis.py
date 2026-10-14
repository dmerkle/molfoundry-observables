"""Discover chlorination counts and certify the exact joint event model."""
from __future__ import annotations
from dataclasses import replace
import sympy as sp
from mfobs.api import ModelInputs, Observable, RatedRule
from mfobs.core.conditions import UnsupportedCondition
from mfobs.core.pattern_support import discover_support, polynomial_record
from mfobs.models.chlorination.inputs import chlorination_rule
from mfobs.models.chlorination.counting import LocalConditionalReducer
from mfobs.models.chlorination.events import _certify_events, _certify_growth, _derive_events


class Chlorination:
    """Drive support and generator images from actual MolFoundry commutators."""

    def __init__(self, variant="guarded", *, chlorine_total=4, state_space=None, rules=None):
        if variant not in ("guarded", "ablated"):
            raise ValueError(f"unknown chlorination variant {variant!r}")
        guarded = variant == "guarded"
        from mfobs.models.chlorination.inputs import state_family as chlorinated_rings
        expected = chlorinated_rings(chlorine_total)
        self.omega = expected if state_space is None else state_space.require(expected)
        if self.omega.definition["chlorine_total"] != 4:
            raise ValueError("the declared family fixes total chlorine at four")
        self.chlorine_total = self.omega.definition["chlorine_total"]
        self.guarded = guarded
        self.rule = chlorination_rule(guarded)
        self.reducer = LocalConditionalReducer(self.chlorine_total, state_space=self.omega)
        self.patterns = self.reducer.patterns
        from mfobs.models.chlorination.inputs import grow_rule, RATE_FACTORS
        self.growth_rule = grow_rule()
        if rules is not None:
            allowed = {chlorination_rule(value).canonKey for value in (True, False)}
            if (set(rules) != {"chlorinate", "grow"} or
                    rules["chlorinate"].canonKey not in allowed or
                    rules["grow"].canonKey != self.growth_rule.canonKey):
                raise ValueError("unsupported rules for the chlorination counting family")
            self.rule, self.growth_rule = rules["chlorinate"], rules["grow"]
            self.guarded = self.rule.matchCondition is not None
        self.omega.validate_rule(self.rule)
        self.omega.validate_rule(self.growth_rule)
        self.growth_rate_factor = RATE_FACTORS["grow"]

    def specification(self) -> ModelInputs:
        return ModelInputs(self.omega,
            {"chlorinate": RatedRule(self.rule, "k", identifier="chlorinate"),
             "grow": RatedRule(self.growth_rule, "g", self.growth_rate_factor)},
            {"Z": Observable("Z", self.patterns["Z"].rule)})

    def rule_gmls(self):
        return self.specification().rule_gmls()

    def query_gmls(self):
        return self.specification().query_gmls()

    def run(self, inputs: ModelInputs, *, maximum_patterns=16,
            maximum_expansions=16):
        inputs.require(self.specification())
        self.rule = inputs.rules["chlorinate"].rule
        self.growth_rule = inputs.rules["grow"].rule
        self.omega = self.reducer.omega = inputs.omega
        _, identity = inputs.observables["Z"].identity_terms()[0]
        self.patterns["Z"] = replace(self.patterns["Z"], rule=identity)
        queries = {"Z": self.patterns["Z"]}
        reaction = inputs.rules["chlorinate"]
        self.reducer.trace = {}
        self.reducer.activity_reduction = None
        support = discover_support(
            (reaction,), queries,
            counter=self.reducer,
            maximum_patterns=maximum_patterns,
            maximum_expansions=maximum_expansions,
        )
        engine_images = {
            name: support.images[name][reaction.identifier]
            for name in support.patterns
            if reaction.identifier in support.images.get(name, {})
        }
        complete = support.status == "POLYNOMIAL_SUPPORT_CLOSED"
        names = tuple(support.patterns)
        events, certificate, event_reason, growth_certificate = (), None, None, None
        if complete:
            try:
                growth_certificate = _certify_growth(self.growth_rule, support.patterns)
                events = _derive_events(self.rule, self.reducer, names)
                certificate = _certify_events(self.rule, events, self.reducer, names)
            except UnsupportedCondition as exc:
                complete = False
                event_reason = str(exc)
                events = ()
        if complete:
            event_activity = sp.expand(sum(e.activity for e in events))
            if not self.reducer.equivalent(event_activity, support.activities[reaction.identifier]):
                raise AssertionError("event activity differs from raw-composition activity")
            for index, name in enumerate(names):
                from_events = sp.expand(sum(e.activity * e.delta[index] for e in events))
                if not self.reducer.equivalent(from_events, engine_images[name]):
                    raise AssertionError(f"event model differs from raw G{name}")

        result = {
            "schema_version": 5,
            "status": "PASS" if complete else "UNKNOWN",
            "rule": {"gml": self.rule.getGMLString(),
                     "condition_kind": "generic-not-exists" if self.guarded else "none"},
            "composition_bridge": {
                "rule_channel": "chlorinate",
                "all_generator_images_from_raw_commutators": True,
                "all_projected_terms_reduced_from_engine_conditions": True,
                "activity_from_dpo_jump_closure": True,
                "state_space_enumeration": False,
                "activity_reduction": self.reducer.activity_reduction.evidence if self.reducer.activity_reduction else None,
                "generator_reductions": self.reducer.trace,
            },
            "synthesis": {
                "selected_query": "Z",
                "counting_vocabulary": "positive patterns constructed from conditional-count monomials",
                "state_language": {
                    "head": "one alkylcyclopentadienyl anion with arbitrary finite linear tail",
                    "sites": "four neutral positions, each bearing H or Cl",
                    "chlorine_total": self.chlorine_total,
                    "HCl_count": "4-Z", "methane_count": "arbitrary finite",
                    "H2_count": "arbitrary finite", "tail_extension": "n >= 0, unbounded",
                    "finite": False,
                    "preservation": "chlorination preserves the tail; growth adds one tail carbon, "
                                    "consumes CH4 and produces H2; both preserve the aromatic backbone and chlorine inventory",
                },
                "state_space_enumeration": False,
                "positive_patterns": {name: p.metadata for name, p in support.patterns.items()},
                "support": support.record(),
            },
            "growth": {
                "rule_gml": self.growth_rule.getGMLString(),
                "per_match_factor": str(self.growth_rate_factor),
                "certificate": growth_certificate,
                "generator_images": {name: "0" for name in names} if growth_certificate else {},
                "projection": "growth is an active molecular channel with zero observable jump; "
                              "its nonzero activity is not set to zero or included in the reduced event table",
            },
            "events": [{
                "id": e.identifier,
                "activity": polynomial_record(e.activity),
                "delta": {n: v for n, v in zip(names, e.delta) if v},
                "certificate": e.certificate,
            } for e in events],
            "pjc": {
                "status": "CLOSED_FINITE_EVENT_MODEL" if complete else "UNKNOWN",
                "coordinates": list(names), "reason": event_reason, "event_count": len(events),
                "all_activities_polynomial": complete, "constant_integer_jumps": complete,
                "generator_match": "exact equality with raw-composition generator images" if complete else None,
                "joint_partition": certificate,
            },
        }
        return result
