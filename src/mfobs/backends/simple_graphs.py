"""Exact positive-pattern backend for finite simple labelled graphs.

This backend is intentionally modest but reusable. It is not a tutorial
special case: any unconditional DPO model over the declared vertex/edge
signature can use it when the required reductions stay inside positive
identity graph counts.

The backend performs no finite-state enumeration. New graph observables are
introduced from actual MolFoundry compositions. Unsupported reductions return
an error rather than silently changing semantics.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from typing import Mapping

import sympy as sp

from mfobs.api import JumpClosureResult, ModelInputs, Observable, RatedRule
from mfobs.core.conditions import engine_has_conditions
from mfobs.core.rules import rule_graph, rule_vertex_ids, added_edges
from mfobs.core.composition_jump_model import (
    GraphReactionClass,
)
from mfobs.core.engine import (
    commutator,
    dpo_jump_closure,
    mf,
)
from mfobs.core.uniform_events import extract_uniform_events, constant_partition, compile_event_model
from mfobs.core.graph import node_rows, edge_rows
from mfobs.core.identity import (
    identity_contradictions,
    jump_input,
    normalized_identity_rule,
)
from mfobs.core.jump_closure_worklist import (
    ObservableClosureResult,
    ObservableDesignation,
    RuleExpression,
    bilinear_commutator,
    discover_observables_by_channel,
    jump_close_expression,
)
from mfobs.core.observables import parse_observable
from mfobs.core.pattern_support import (
    Image,
    Pattern,
    UnsupportedReduction,
    discover_support,
)
from mfobs.core.polynomial_jump_closure import (
    PolynomialRing,
    JumpModelProvenance,
    analyze_polynomial_jump_closure,
)


@dataclass(frozen=True)
class _Side:
    raw: int
    kept: int

    def record(self) -> dict:
        return {"raw": self.raw, "kept": self.kept}


@dataclass(frozen=True)
class _RawCommutator:
    coefficients: Mapping[str, int]
    representatives: Mapping[str, object]
    left: _Side
    right: _Side


class _PatternRegistry:
    """Stable readable names for identity graphs discovered at runtime."""

    def __init__(self):
        self.key_to_name: dict[str, str] = {}
        self.name_to_key: dict[str, str] = {}
        self.counter = 0

    @staticmethod
    def _suggest(rule) -> str | None:
        parsed = parse_observable(normalized_identity_rule(rule))
        if parsed.no_edges or parsed.generic_conditions:
            return None
        labels = parsed.labels
        if len(labels) == 1 and not parsed.edges:
            return "N_" + next(iter(labels.values()))
        if len(labels) == 2 and len(parsed.edges) == 1:
            u, v, _ = parsed.edges[0]
            a, b = sorted((labels[u], labels[v]))
            return "E_" + a + b
        return None

    def register(self, rule, *, preferred: str | None = None) -> str:
        normalized = normalized_identity_rule(rule)
        key = normalized.canonKey
        if key in self.key_to_name:
            return self.key_to_name[key]

        candidate = preferred or self._suggest(normalized)
        if candidate is None or (
            candidate in self.name_to_key and self.name_to_key[candidate] != key
        ):
            while True:
                candidate = f"O{self.counter}"
                self.counter += 1
                if candidate not in self.name_to_key:
                    break

        self.key_to_name[key] = candidate
        self.name_to_key[candidate] = key
        return candidate

    def make_pattern(self, name: str, rule, kind: str) -> Pattern:
        normalized = SimpleGraphBackend._positive_identity(rule)
        return Pattern(
            name=name,
            rule=normalized,
            kind=kind,
            metadata={"canonical_key": normalized.canonKey},
        )

    def activity(self, reaction: GraphReactionClass) -> Image:
        identity = SimpleGraphBackend._positive_identity(
            dpo_jump_closure(
                reaction.rule,
                name=f"activity_{reaction.identifier}",
            )
        )
        name = self.register(identity)
        return Image(
            expression=sp.Symbol(name),
            factors={name: self.make_pattern(name, identity, "activity-graph")},
            evidence={
                "source": "DPO activity identity of the literal graph rule",
                "activity_identity_gml": identity.getGMLString(),
            },
        )

    def image(
        self, pattern: Pattern,
        reaction: GraphReactionClass,
    ) -> Image:
        raw, pair_records = bilinear_commutator(
            RuleExpression.from_rules((pattern.rule,)),
            RuleExpression.from_rules((reaction.rule,)),
            SimpleGraphBackend._pair,
        )
        jumped, jump_trace = jump_close_expression(
            raw,
            jump_closure=dpo_jump_closure,
        )

        expression = sp.Integer(0)
        factors = {}
        for key, coefficient in jumped.coefficients.items():
            identity = SimpleGraphBackend._positive_identity(
                jumped.representatives[key]
            )
            if identity_contradictions(identity):
                continue
            name = self.register(identity)
            value = Fraction(coefficient)
            expression += sp.Rational(
                value.numerator,
                value.denominator,
            ) * sp.Symbol(name)
            factors[name] = self.make_pattern(
                name,
                identity,
                "generator-image-graph",
            )

        return Image(
            expression=sp.expand(expression),
            factors=factors,
            evidence={
                "source": "MolFoundry composition plus DPO jump projection",
                "raw_commutator": raw.record(),
                "composition_pairs": pair_records,
                "jump_closure_trace": jump_trace,
            },
        )



class SimpleGraphBackend:
    """Backend for all finite simple undirected graphs over a fixed signature.

    It supports unconditional DPO rules and positive identity-pattern counts.
    Omega consists of all finite simple graphs using only the declared vertex
    and edge labels. The rule-closure check is conservative: an edge created
    between two preserved vertices is rejected because, without a NoEdge
    condition, some host in Omega could already contain that edge.

    Algorithm 2 additionally needs a uniform-match argument. The public
    jump_closure call therefore requires uniform_matches=True before event
    extraction. This is an explicit semantic assertion, not inferred from
    aggregate generator images.
    """

    def __init__(self, omega):
        from mfobs.state_space import simple_graphs
        omega.require(simple_graphs(vertex_labels=omega.definition["vertex_labels"],
                                    edge_labels=omega.definition["edge_labels"]))
        self.omega = omega
        self.vertex_labels = frozenset(omega.definition["vertex_labels"])
        self.edge_labels = frozenset(omega.definition["edge_labels"])

    def _contains(self, graph: mf.Graph) -> bool:
        return (
            all(label in self.vertex_labels for _, label in node_rows(graph))
            and all(label in self.edge_labels for _, _, label in edge_rows(graph))
        )

    def _check_rule_closed(self, declaration: RatedRule) -> None:
        if declaration.per_match != 1:
            raise ValueError("SimpleGraphs requires unit per-match factors")
        rule = declaration.rule
        if not self._contains(rule_graph(rule, 'left')) or not self._contains(rule_graph(rule, 'right')):
            raise ValueError(f"{rule.name}: rule uses labels outside Omega")
        if rule_vertex_ids(rule, 'changed'):
            raise ValueError(f"{rule.name}: relabelling is not supported by SimpleGraphs")

        fresh = set(rule_vertex_ids(rule, 'created'))
        for u, v, _ in added_edges(rule):
            if u not in fresh and v not in fresh:
                raise ValueError(
                    f"{rule.name}: created edge {u}-{v} joins two preserved "
                    "vertices; SimpleGraphs cannot prove rule-closure on all "
                    "simple hosts without an application condition"
                )
        if engine_has_conditions(rule):
            raise ValueError(f"{rule.name}: SimpleGraphs supports unconditional rules only")

    @staticmethod
    def _engine_observable(observable: Observable):
        return mf.Rule.identity(observable.pattern, observable.name, add=False)

    @staticmethod
    def _pair(observable, reaction) -> _RawCommutator:
        result = commutator(
            observable,
            reaction,
            connected=False,
            include_empty=True,
        )
        return _RawCommutator(
            coefficients=result.composite_coefficients,
            representatives=result.composite_representatives,
            left=_Side(result.left_raw, result.left_raw),
            right=_Side(result.right_raw, result.right_raw),
        )

    @staticmethod
    def _positive_identity(rule):
        _, deleted = jump_input(rule)
        if deleted:
            # On all simple graphs a deleted vertex may carry further edges, so
            # the DPO dangling condition cannot be dropped.
            raise UnsupportedReduction(
                f"DPO jump closure deletes vertices {list(deleted)}: its dangling "
                "condition is not a positive pattern count"
            )
        normalized = normalized_identity_rule(rule)
        parsed = parse_observable(normalized)
        if parsed.no_edges or parsed.generic_conditions:
            raise UnsupportedReduction(
                "SimpleGraphs currently supports positive identity patterns only"
            )
        return normalized

    def _prepare(self, inputs: ModelInputs):
        """Validate native inputs and start a query-only pattern registry."""
        inputs.omega.require(self.omega)
        if len(inputs.observables) != 1:
            raise ValueError("SimpleGraphBackend requires one selected graph query")
        rules = tuple(inputs.rules.values())
        query = next(iter(inputs.observables.values()))
        if not isinstance(query.pattern, mf.Graph) or query.terms:
            raise ValueError("SimpleGraphs requires a single positive graph query")
        if not self._contains(query.pattern):
            raise ValueError("query pattern lies outside Omega")
        for rule in rules:
            self._check_rule_closed(rule)

        seed_rule = self._engine_observable(query)
        registry = _PatternRegistry()
        registry.register(seed_rule, preferred=query.name)

        return rules, seed_rule, registry

    def linear_closure(
        self,
        inputs: ModelInputs,
        *,
        maximum_dimension: int = 16,
    ) -> ObservableClosureResult:
        """Algorithm 1: discover an invariant linear span, by rate channel."""
        rules, seed_rule, registry = self._prepare(inputs)
        query = next(iter(inputs.observables.values()))
        channel_rules: defaultdict[str, list[object]] = defaultdict(list)
        for declaration in rules:
            channel_rules[declaration.rate].append(declaration.rule)

        generators = {
            channel: RuleExpression.from_rules(group)
            for channel, group in channel_rules.items()
        }

        def classify(expression: RuleExpression, index: int):
            if len(expression.coefficients) == 1:
                key, coefficient = next(iter(expression.coefficients.items()))
                if coefficient == 1:
                    rule = expression.representatives[key]
                    try:
                        name = registry.register(rule)
                    except ValueError:
                        return ObservableDesignation(
                            f"L{index}",
                            "discovered DPO jump closure with a dangling condition",
                        )
                    return ObservableDesignation(
                        name,
                        "count of the discovered identity graph",
                    )
            return ObservableDesignation(
                f"L{index}",
                "discovered linear combination of identity graph counts",
            )

        return discover_observables_by_channel(
            seeds={query.name: RuleExpression.from_rules((seed_rule,))},
            generators=generators,
            pair_commutator=self._pair,
            classify=classify,
            jump_closure=dpo_jump_closure,
            jump_projection_name="DPO",
            maximum_dimension=maximum_dimension,
        )

    def jump_closure(
        self,
        inputs: ModelInputs,
        *,
        uniform_matches: bool = False,
        maximum_support_patterns: int = 16,
        maximum_support_expansions: int = 16,
    ) -> JumpClosureResult:
        """Algorithm 2: discover pattern support and extract polynomial events.

        ``uniform_matches=True`` asserts the caller's all-state uniform-match
        argument. Without it, support discovery runs but no events are certified.
        """
        rules, seed_rule, registry = self._prepare(inputs)
        query = next(iter(inputs.observables.values()))
        reaction_classes = tuple(
            GraphReactionClass(
                identifier=declaration.identifier,
                channel=declaration.rate,
                rule=declaration.rule,
                source_rule=declaration.identifier,
                predicate=f"one injective match of {declaration.identifier}",
                metadata={"backend": "SimpleGraphs"},
            )
            for declaration in rules
        )

        support = discover_support(
            rules=reaction_classes,
            seeds={
                query.name: registry.make_pattern(
                    query.name,
                    seed_rule,
                    "selected-query",
                )
            },
            counter=registry,
            maximum_patterns=maximum_support_patterns,
            maximum_expansions=maximum_support_expansions,
        )
        if support.status != "POLYNOMIAL_SUPPORT_CLOSED" or not uniform_matches:
            status = (support.status if support.status != "POLYNOMIAL_SUPPORT_CLOSED"
                      else "UNKNOWN_UNIFORMITY")
            return JumpClosureResult(
                support, None,
                {"status": status, "reason": (
                    "Support remains open" if status != "UNKNOWN_UNIFORMITY" else
                    "Event extraction needs a supplied all-state uniform-match justification")},
            )

        # One constant-jump class per source rule, subject to the caller's
        # uniform-match argument. Nonconstant quotients remain unsupported.
        from types import SimpleNamespace
        sources = tuple(SimpleNamespace(identifier=r.identifier, rate=r.channel)
                        for r in reaction_classes)
        argument = "caller-supplied all-state uniform-match justification for each source rule"
        extracted = extract_uniform_events(support, sources,
            partition=constant_partition, scope="declared finite simple labelled graphs",
            uniform_match_argument=argument,
            partition_argument="one exhaustive enabled class per rule; quotient must be constant")
        if extracted["status"] != "FINITE_EVENTS":
            return JumpClosureResult(support, None,
                {"status": extracted["status"], "reason": extracted.get("reason", "")})
        ring = PolynomialRing(tuple(support.patterns))
        digest = sha256(
            (
                json.dumps(inputs.omega.record(), sort_keys=True)
                + query.name
                + query.pattern.getGMLString()
                + "".join(
                    declaration.identifier
                    + declaration.rate
                    + declaration.rule.getGMLString()
                    for declaration in sorted(rules, key=lambda rule: rule.identifier)
                )
            ).encode()
        ).hexdigest()

        compilation = compile_event_model(extracted,
            name="simple-graph-analysis", ring=ring,
            certificate={"status": "VERIFIED_ALL_STATES",
                "argument_status": "CALLER_ASSERTED_UNIFORM_MATCHES",
                "uniform_match_obligations": extracted["uniform_match_obligations"],
                "partition": extracted["partition_argument"],
                "enumeration_used": False},
            input_definition_sha256=digest,
            provenance=JumpModelProvenance(
                backend="SimpleGraphs positive identity-count reduction",
                event_partition_origin=extracted["partition_argument"],
                jump_origin="composition-derived generator images divided by activity",
                activity_origin="DPO activity identity of the supplied source rule",
                jump_projection="DPO", composition_scope="all discovered coordinates and source rules",
                all_jumps_composition_derived=True, supplied_delta_table=False))

        return JumpClosureResult(
            support, compilation, analyze_polynomial_jump_closure(compilation.model),
        )

    def run(
        self,
        inputs: ModelInputs,
        *,
        uniform_matches: bool = False,
        maximum_dimension: int = 16,
        maximum_support_patterns: int = 16,
        maximum_support_expansions: int = 16,
    ) -> dict:
        """Run the two algorithms independently from the same inputs."""
        linear = self.linear_closure(
            inputs, maximum_dimension=maximum_dimension,
        )
        jump = self.jump_closure(
            inputs, uniform_matches=uniform_matches,
            maximum_support_patterns=maximum_support_patterns,
            maximum_support_expansions=maximum_support_expansions,
        )
        events = [] if jump.compilation is None else [
            {"identifier": event.identifier, "rate": event.channel,
             "jump": event.delta, "activity": jump.compilation.model.ring.record(event.propensity)}
            for event in jump.compilation.model.events
        ]
        return {"linear": linear.record(), "support_discovery": jump.support.record(),
                "pjc": jump.pjc, "events": events}



__all__ = ["SimpleGraphBackend"]
