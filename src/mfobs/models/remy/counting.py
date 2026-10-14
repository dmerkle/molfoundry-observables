"""Tree count reduction and linear closure from the selected spine query."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
import json
from mfobs.core.conditions import engine_has_conditions
from mfobs.core.engine import dpo_jump_closure, mf
from mfobs.core.jump_closure_worklist import (
    ObservableDesignation,
    RuleExpression,
    discover_observables,
)
from mfobs.core.rules import rule_graph, rule_vertex_ids, removed_edges, preserves_vertices
from mfobs.core.graph import node_labels, edge_rows
from mfobs.core.typed_incidence import pair_commutator
from mfobs.models.remy.inputs import EDGE_TYPES, TYPED_INCIDENCE_SPECIFICATION, _gml_body, insertion_rules, observable_gml


REMY_JUMP_PROJECTION = "DPO with valid-tree input-count reduction"


def valid_tree_dpo_jump_closure(rule, name: str = "jump"):
    """Reduce the DPO projection's count to a full-input identity on Omega.

    First construct I <- K -> I, retaining DPO deletion semantics. Only
    complete typed edge-nodes may be deleted, with their src/tgt incidences.
    Such nodes have no other incidences in a valid tree, so every injective
    input match satisfies dangling. Recreating them returns an isomorphic
    host. Thus this projection and the full-input identity have equal counts
    on Omega, although their spans differ. Identity-query compositions obey
    the same invariant; it is checked again for every projected term.

    This narrow reducer fails closed outside the unconditioned incidence
    fragment. It is not a generic DPO-to-identity or DPO-to-SqPO conversion.
    """
    projected = dpo_jump_closure(rule, name)
    if engine_has_conditions(projected):
        raise ValueError("Rémy count reduction does not support additional conditions")
    # The projection deletes exactly what the rule deletes.  Read that from the
    # rule: an edge recreated between kept vertices is kept in the projection.
    left = rule_graph(rule)
    nodes, edges = node_labels(left), edge_rows(left)
    if set(nodes.values()) - {"V", *EDGE_TYPES}:
        raise ValueError("Rémy count reduction requires typed-incidence nodes")
    deleted = set(rule_vertex_ids(rule, "deleted"))
    for node in deleted:
        incident = [(u, v, label) for u, v, label in edges if node in (u, v)]
        if (nodes[node] not in EDGE_TYPES
                or sorted(label for _, _, label in incident) != ["src", "tgt"]
                or any(nodes[v if u == node else u] != "V" for u, v, _ in incident)):
            raise ValueError("DPO dangling is not certified: deleted node needs complete src/tgt incidences")
    if any(u not in deleted and v not in deleted for u, v, _ in removed_edges(rule)):
        raise ValueError("Rémy count reduction only supports deletion of complete edge-nodes")
    gml = ('rule [ ruleID ' + json.dumps(name) + '\n context [\n'
           + _gml_body(nodes.items(), edges) + '\n ]\n]')
    return mf.Rule.fromGMLString(gml, add=False)


def _left_spine_component(rule) -> tuple[int, str]:
    """Recognise one identity component by its directed incidence structure.

    This is a post-discovery labeller, not a source of candidate observables.
    It reconstructs the represented directed typed graph, identifies its
    exterior stem, and verifies the full ordered-binary left-spine shape.
    """

    if not preserves_vertices(rule):
        raise ValueError("an observable component must be an identity rule")
    kept = rule_graph(rule, "context")
    labels = node_labels(kept)
    incidence = {
        node: {"src": [], "tgt": []}
        for node, label in labels.items()
        if label in EDGE_TYPES
    }
    for first, second, bond in edge_rows(kept):
        edge_nodes = [node for node in (first, second) if node in incidence]
        if len(edge_nodes) != 1 or bond not in {"src", "tgt"}:
            raise ValueError("not a typed-incidence identity pattern")
        edge_node = edge_nodes[0]
        object_node = second if first == edge_node else first
        if labels.get(object_node) != "V":
            raise ValueError("typed edge is not incident to an object node")
        incidence[edge_node][bond].append(object_node)
    represented = []
    for edge_node, ports in incidence.items():
        if len(ports["src"]) != 1 or len(ports["tgt"]) != 1:
            raise ValueError("typed edge must have exactly one source and target")
        represented.append(
            (ports["src"][0], ports["tgt"][0], labels[edge_node], edge_node)
        )

    objects = {node for node, label in labels.items() if label == "V"}
    outgoing = {source: edge for source, _, _, edge in represented}
    exteriors = objects - set(outgoing)
    if len(exteriors) != 1:
        raise ValueError("left-spine component must have one exterior object")
    exterior = next(iter(exteriors))
    stems = [edge for edge in represented if edge[1] == exterior]
    if len(stems) != 1:
        raise ValueError("left-spine component must have one exterior stem")
    current, _, stem_type, stem_node = stems[0]
    if len(represented) == 1:
        if objects != {current, exterior}:
            raise ValueError("edge observable has unexpected objects")
        return 0, stem_type

    incoming: dict[int, dict[str, int]] = defaultdict(dict)
    for source, target, edge_type, edge_node in represented:
        if edge_node == stem_node:
            continue
        if edge_type not in {"L", "R"} or edge_type in incoming[target]:
            raise ValueError("not an ordered full-binary pattern")
        incoming[target][edge_type] = source
    visited = {exterior}
    depth = 0
    while True:
        children = incoming.get(current, {})
        if set(children) != {"L", "R"}:
            raise ValueError("left spine ended at a non-full internal vertex")
        visited.add(current)
        right = children["R"]
        if right in incoming:
            raise ValueError("right branch is not a leaf in a left-spine pattern")
        visited.add(right)
        depth += 1
        left = children["L"]
        if left not in incoming:
            visited.add(left)
            break
        current = left
    if visited != objects:
        raise ValueError("left-spine component has additional structure")
    return depth, stem_type


def classify_observable(
    expression: RuleExpression, index: int
) -> ObservableDesignation:
    """Assign the conventional P_k name after algebraic discovery."""

    if set(expression.coefficients.values()) != {Fraction(1)}:
        raise ValueError("Rémy observable components must have unit weight")
    shapes = [
        _left_spine_component(expression.representatives[key])
        for key in expression.coefficients
    ]
    depths = {depth for depth, _ in shapes}
    stems = {stem for _, stem in shapes}
    if len(shapes) != 3 or len(depths) != 1 or stems != set(EDGE_TYPES):
        raise ValueError("computed residual is not a complete typed P_k family")
    depth = next(iter(depths))
    name = "E" if depth == 0 else f"P{depth}"
    meaning = (
        "all selectable directed typed edges"
        if depth == 0
        else f"full ordered-binary left-spine pattern of depth {depth}"
    )
    return ObservableDesignation(
        name,
        meaning,
        {
            "classification": "post-discovery structural recognition",
            "left_spine_depth": depth,
            "stem_components": sorted(stems),
            "generic_fallback_name": f"O{index}",
        },
    )


def observable_synthesis(order: int = 3, *, maximum_dimension: int = 64):
    """Discover the first-commutator observable span from one P_k alone.

    The generic core creates every residual expression.  This model adapter
    supplies the six rules, the declared incidence/acyclicity theory, and a
    post-hoc structural name for each expression; it supplies no target basis
    and no commutator equation.
    """

    rules = insertion_rules()
    seed_name = f"P{order}"
    seed_components = tuple(
        mf.Rule.fromGMLString(observable_gml(seed_name, stem), add=False)
        for stem in EDGE_TYPES
    )
    seed = RuleExpression.from_rules(seed_components)
    generator = RuleExpression.from_rules(rules.values())

    result = discover_observables(
        seeds={seed_name: seed},
        generator=generator,
        pair_commutator=pair_commutator(TYPED_INCIDENCE_SPECIFICATION, include_empty=True),
        classify=classify_observable,
        jump_closure=valid_tree_dpo_jump_closure,
        jump_projection_name=REMY_JUMP_PROJECTION,
        maximum_dimension=maximum_dimension,
    )
    result.origins[seed_name].update(
        meaning=f"full ordered-binary left-spine pattern of depth {order}",
        metadata={
            "left_spine_depth": order,
            "stem_components": list(EDGE_TYPES),
        },
    )
    return result
