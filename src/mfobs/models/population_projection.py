"""Generic typed-token projection for finite-component population rules.

An application supplies molecular component inventories and classifies each
component name as dynamic, persistent, or discarded.  This module constructs
small typed DPO rules mechanically.  Observable jumps are deliberately absent:
they are recovered later by the common composition jump compiler.
"""
from __future__ import annotations

from mfobs.core.graph import node_rows, make_graph, edge_rows, node_labels
from mfobs.core.rules import make_rule
from mfobs.core.engine import mf
from collections import Counter
from dataclasses import dataclass
from typing import Mapping

from mfobs.core.composition_jump_model import GraphReactionClass
from mfobs.core.gml import component_graphs

from mfobs.core.jump_closure_worklist import RuleExpression
from mfobs.core.molfoundry_execution import canonical_key
from mfobs.core.observables import parse_observable
from mfobs.core.polynomial_jump_closure import Polynomial, PolynomialRing
from mfobs.models.aba.projection import KIND_BOND, POSITION, TerminalWordGraphModel, TypedGraphBuilder


def molecular_component_inventory(
    graph: mf.Graph,
    catalogue: Mapping[str, mf.Graph],
) -> dict[str, int]:
    """Read a complete-component inventory using a declared finite catalogue."""

    keys = {canonical_key(value): name for name, value in catalogue.items()}
    if len(keys) != len(catalogue):
        raise ValueError("component catalogue entries must be pairwise distinct")
    result = Counter()
    for component in component_graphs(graph):
        key = canonical_key(component)
        if key not in keys:
            raise ValueError("molecular rule side contains an undeclared component")
        result[keys[key]] += 1
    return dict(sorted(result.items()))


def _falling_factorial(
    ring: PolynomialRing,
    coordinate: str,
    order: int,
) -> Polynomial:
    value = ring.constant(1)
    variable = ring.coordinate(coordinate)
    for offset in range(order):
        value = ring.multiply(
            value,
            ring.subtract(variable, ring.constant(offset)),
        )
    return value


@dataclass(frozen=True)
class PopulationChannel:
    """One molecular channel after its component inventories were read."""

    identifier: str
    channel: str
    left: Mapping[str, int]
    right: Mapping[str, int]


class TypedTokenGraphModel:
    """Mechanically project component stoichiometry to typed token rules."""

    def __init__(
        self,
        schema: TerminalWordGraphModel,
        *,
        dynamic_coordinates: Mapping[str, str],
        persistent_components: tuple[str, ...] = (),
        discarded_components: tuple[str, ...] = (),
    ) -> None:
        self.schema = schema
        self.dynamic_coordinates = dict(dynamic_coordinates)
        self.persistent_components = frozenset(persistent_components)
        self.discarded_components = frozenset(discarded_components)
        groups = (
            set(self.dynamic_coordinates),
            set(self.persistent_components),
            set(self.discarded_components),
        )
        if any(groups[index] & groups[other]
               for index in range(3) for other in range(index)):
            raise ValueError("component roles must be disjoint")

    @staticmethod
    def tag(component: str) -> str:
        return f"Population:{component}"

    @property
    def tags(self) -> tuple[str, ...]:
        names = (
            set(self.dynamic_coordinates)
            | set(self.persistent_components)
            | set(self.discarded_components)
        )
        return tuple(self.tag(name) for name in sorted(names))

    def _add_token(self, builder: TypedGraphBuilder, component: str) -> None:
        self.schema.tagged_position(builder, self.tag(component))

    def token_graph(self, component: str, count: int = 1) -> mf.Graph:
        """Build a disjoint typed-token graph for reusable rule context."""

        if component not in self.dynamic_coordinates:
            raise ValueError("token context must name a dynamic component")
        if count < 0:
            raise ValueError("token count must be nonnegative")
        builder = TypedGraphBuilder.empty()
        for _ in range(count):
            self._add_token(builder, component)
        return builder.graph()

    def structural_rule(self, channel: PopulationChannel) -> mf.Rule:
        """Build the DPO rule implied by one classified component inventory."""

        known = (
            set(self.dynamic_coordinates)
            | set(self.persistent_components)
            | set(self.discarded_components)
        )
        if (set(channel.left) | set(channel.right)) - known:
            raise ValueError(f"{channel.identifier}: unclassified component")

        common = TypedGraphBuilder.empty()
        left_only = Counter()
        right_only = Counter()
        for component in sorted(known):
            before = int(channel.left.get(component, 0))
            after = int(channel.right.get(component, 0))
            if min(before, after) < 0:
                raise ValueError("component inventories must be nonnegative")
            if component in self.persistent_components:
                if before > 1 or after > 1:
                    raise ValueError("persistent components have unit population")
                if before or after:
                    self._add_token(common, component)
                continue
            if component in self.discarded_components:
                if before:
                    raise ValueError("discarded components may occur only as products")
                continue
            retained = min(before, after)
            for _ in range(retained):
                self._add_token(common, component)
            left_only[component] = before - retained
            right_only[component] = after - retained

        left = self.schema.copy_graph(common.graph())
        for component, count in sorted(left_only.items()):
            for _ in range(count):
                self._add_token(left, component)
        right = self.schema.copy_graph(common.graph())
        right.next_id = max(left.nodes, default=-1) + 1
        for component, count in sorted(right_only.items()):
            for _ in range(count):
                self._add_token(right, component)
        return make_rule(channel.identifier, left.graph(), right.graph())

    def reaction_class(
        self, channel: PopulationChannel
    ) -> GraphReactionClass:
        return GraphReactionClass(
            channel.identifier,
            channel.channel,
            self.structural_rule(channel),
            channel.identifier,
            "component inventory class derived from the molecular rule sides",
            {
                "kind": "finite-component-population",
                "molecular_left_inventory": dict(channel.left),
                "molecular_right_inventory": dict(channel.right),
            },
        )

    def observable_expression(self, component: str) -> RuleExpression:
        if component not in self.dynamic_coordinates:
            raise ValueError("only dynamic components have population observables")
        builder = TypedGraphBuilder.empty()
        self._add_token(builder, component)
        identity = make_rule(component, builder.graph(), builder.graph())
        return RuleExpression.from_rules((
            identity,
        ))

    def dangling_certificate(self, fallback=None):
        """Certify DPO dangling of deleted token positions; delegate the rest.

        Declared language invariant: a token is exactly one position with one
        kind relation to its population tag.  Token rules create and delete
        tokens only as whole isolated components.  Word rules attach
        next/terminal incidences only to the chain end and to a fed position
        whose kind tag is a word symbol; a position has one kind relation and
        only such fed positions become chain ends, so neither is ever a token.
        Hence a deleted position whose kind relation to a population tag lies
        in ``I`` has no other host edge.  The activity reduction below relies
        on the same invariant.
        """
        from mfobs.core.identity import DanglingNotCertified

        tags = set(self.tags)

        def certify(graph: mf.Graph, deleted) -> None:
            labels = node_labels(graph)
            rest = []
            for vertex in deleted:
                kinds = [labels[v if u == vertex else u]
                         for u, v, bond in edge_rows(graph)
                         if vertex in (u, v) and bond == KIND_BOND]
                if labels[vertex] == POSITION and any(tag in tags for tag in kinds):
                    continue
                rest.append(vertex)
            if rest:
                if fallback is None:
                    raise DanglingNotCertified(
                        f"deleted vertices {rest} are not token positions")
                fallback(graph, tuple(rest))

        return certify

    def activity_polynomial(
        self,
        ring: PolynomialRing,
        _reaction: GraphReactionClass,
        activity_rule,
    ) -> Polynomial:
        """Reduce a generated token activity identity by factorial counting."""

        parsed = parse_observable(activity_rule)
        if parsed.no_edges or parsed.generic_conditions:
            raise ValueError("token activities must be unconditional identities")
        graph = make_graph(parsed.labels, parsed.edges)
        inventory = Counter()
        for component in component_graphs(graph):
            labels = Counter(node_labels(component).values())
            if labels[POSITION] != 1 or len(node_rows(component)) != 2:
                raise ValueError("activity component is not one typed token")
            edges = tuple(edge_rows(component))
            if len(edges) != 1 or edges[0][2] != KIND_BOND:
                raise ValueError("typed token lacks its kind relation")
            tags = [label for label in labels if label != POSITION]
            if len(tags) != 1 or not tags[0].startswith("Population:"):
                raise ValueError("unknown token activity tag")
            inventory[tags[0].split(":", 1)[1]] += 1

        factors = []
        for component, count in sorted(inventory.items()):
            if component in self.dynamic_coordinates:
                factors.append(_falling_factorial(
                    ring, self.dynamic_coordinates[component], count
                ))
            elif component in self.persistent_components:
                if count != 1:
                    return ring.constant(0)
            else:
                raise ValueError("discarded token occurred in an activity identity")
        return ring.multiply(*factors)

    def projection_record(self, channels: tuple[PopulationChannel, ...]) -> dict:
        return {
            "schema_version": 1,
            "kind": "typed-finite-component-population-DPO-projection",
            "dynamic_coordinates": dict(self.dynamic_coordinates),
            "persistent_components": sorted(self.persistent_components),
            "discarded_components": sorted(self.discarded_components),
            "channels": [
                {
                    "identifier": value.identifier,
                    "left_inventory": dict(value.left),
                    "right_inventory": dict(value.right),
                }
                for value in channels
            ],
            "construction": (
                "Molecular side inventories are read by complete-component "
                "isomorphism and translated mechanically to typed DPO tokens; "
                "no observable jump vector is an input."
            ),
        }
