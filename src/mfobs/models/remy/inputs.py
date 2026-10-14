"""Rémy model: rooted planar binary trees, six incidence rules and selected spine counts."""
from __future__ import annotations
from fractions import Fraction
from typing import Iterable
from mfobs.core.engine import mf
from mfobs.state_space import StateSpace


def state_family():
    """Declare Ω for this model, independently of its selected queries."""
    incidence = {
        "name": "directed-planar-binary-forest-negative-patterns-v1",
        "objectLabels": ["V"], "edgeNodeLabels": ["I", "L", "R"],
        "sourceBond": "src", "targetBond": "tgt", "atMostOneOutgoing": True,
        "atMostOneIncomingTypes": ["I", "L", "R"],
        "forbidSuccessorAfterTypes": ["I"], "simpleDirected": True,
        "simpleUnderlying": True, "acyclicDirected": True,
    }
    return StateSpace("planar-binary-trees", {
        "trees": "one finite rooted full planar binary tree",
        "root": "one I stem from the tree root to an exterior vertex",
        "children": ["L", "R"], "orientation": "child edges point toward the tree root",
        "incidence": incidence,
    })


Tree = tuple["Tree", "Tree"] | tuple


LEAF: Tree = ()


EDGE_TYPES = ("I", "L", "R")


OBSERVABLES = ("E", "P1", "P2", "P3")


VARIABLES = ("epsilon", "gamma", "mu", "nu")


def observable_names(order: int = 3) -> tuple[str, ...]:
    if order < 1:
        raise ValueError("Rémy analysis requires a positive P_k order")
    return ("E", *(f"P{depth}" for depth in range(1, order + 1)))


OMEGA = state_family()


TYPED_INCIDENCE_SPECIFICATION = OMEGA.record()["definition"]["incidence"]


def _gml_body(nodes: Iterable[tuple[int, str]], edges: Iterable[tuple[int, int, str]]) -> str:
    lines = [
        f'    node [ id {node} label "{label}" ]'
        for node, label in nodes
    ]
    lines.extend(
        f'    edge [ source {source} target {target} label "{label}" ]'
        for source, target, label in edges
    )
    return "\n".join(lines)


def insertion_rule_gml(side: str, stem_type: str) -> str:
    """Return one of the six rules ``G_side,stem_type``.

    A typed edge is an edge-node with ``src`` and ``tgt`` incidence bonds.
    The original edge-node is deleted and three new typed edge-nodes are
    created; both original object endpoints are preserved.
    """

    if side not in {"L", "R"} or stem_type not in EDGE_TYPES:
        raise ValueError("side must be L/R and stem type must be I/L/R")
    context_nodes = [(0, "V"), (1, "V")]
    left_nodes = [(2, stem_type)]
    left_edges = [(0, 2, "src"), (2, 1, "tgt")]
    right_nodes = [(3, "V"), (4, "V"), (5, stem_type), (6, "L"), (7, "R")]
    if side == "L":
        branches = [
            (0, 6, "src"), (6, 3, "tgt"),
            (4, 7, "src"), (7, 3, "tgt"),
        ]
    else:
        branches = [
            (4, 6, "src"), (6, 3, "tgt"),
            (0, 7, "src"), (7, 3, "tgt"),
        ]
    right_edges = [(3, 5, "src"), (5, 1, "tgt"), *branches]
    return f'''rule [
  ruleID "G{side}_{stem_type}"
  left [
{_gml_body(left_nodes, left_edges)}
  ]
  context [
{_gml_body(context_nodes, [])}
  ]
  right [
{_gml_body(right_nodes, right_edges)}
  ]
]
'''


def insertion_rules() -> dict[str, object]:
    return {
        f"G{side}_{stem}": mf.Rule.fromGMLString(
            insertion_rule_gml(side, stem), add=False
        )
        for side in ("L", "R")
        for stem in EDGE_TYPES
    }


def _observable_typed_edges(name: str, stem_type: str) -> tuple[range, list[tuple[int, int, str]]]:
    if name == "E":
        return range(2), [(0, 1, stem_type)]
    if not name.startswith("P") or not name[1:].isdigit() or int(name[1:]) < 1:
        raise ValueError(f"unknown observable: {name}")
    depth = int(name[1:])
    # A P_k occurrence is an internal root with a left spine of k internal
    # vertices.  Every internal vertex includes both ordered children; the
    # bottom observable root also includes its outgoing edge of type T.
    edges: list[tuple[int, int, str]] = []
    next_vertex = 0
    child = next_vertex
    sibling = next_vertex + 1
    parent = next_vertex + 2
    edges.extend([(child, parent, "L"), (sibling, parent, "R")])
    next_vertex = 3
    current = parent
    for _ in range(1, depth):
        sibling = next_vertex
        parent = next_vertex + 1
        edges.extend([(current, parent, "L"), (sibling, parent, "R")])
        next_vertex += 2
        current = parent
    outside = next_vertex
    edges.append((current, outside, stem_type))
    return range(outside + 1), edges


def observable_gml(name: str, stem_type: str) -> str:
    """Identity-rule component of E/P1/P2/P3 for outgoing type ``stem``."""

    vertices, typed_edges = _observable_typed_edges(name, stem_type)
    nodes = [(vertex, "V") for vertex in vertices]
    incidence_edges = []
    edge_id = len(nodes)
    for source, target, edge_type in typed_edges:
        nodes.append((edge_id, edge_type))
        incidence_edges.extend(
            [(source, edge_id, "src"), (edge_id, target, "tgt")]
        )
        edge_id += 1
    return f'''rule [
  ruleID "{name}_{stem_type}"
  context [
{_gml_body(nodes, incidence_edges)}
  ]
]
'''


def observable_rules(order: int = 3) -> dict[str, tuple[object, ...]]:
    return {
        name: tuple(
            mf.Rule.fromGMLString(observable_gml(name, stem), add=False)
            for stem in EDGE_TYPES
        )
        for name in observable_names(order)
    }


def query_order(selection) -> int:
    """Validate one common file-based query as a complete typed P_k sum."""

    from mfobs.models.remy.counting import _left_spine_component

    if len(selection.observables) != 1:
        raise ValueError("Rémy accepts exactly one selected P_k expression")
    query = selection.observables[0]
    if any(term.coefficient != Fraction(1) for term in query.terms):
        raise ValueError("Rémy's three stem components must have coefficient 1")
    shapes = [
        _left_spine_component(mf.Rule.fromGMLString(term.gml, add=False))
        for term in query.terms
    ]
    depths = {depth for depth, _ in shapes}
    stems = {stem for _, stem in shapes}
    if len(shapes) != 3 or len(depths) != 1 or stems != set(EDGE_TYPES):
        raise ValueError(
            "Rémy P_k must be the sum of exactly the I, L and R stem components"
        )
    order = next(iter(depths))
    if order < 1 or query.name != f"P{order}":
        raise ValueError(f"canonical name for this Rémy expression is P{order}")
    return order


def _typed_tree_graph(tree: Tree) -> tuple[list[tuple[int, str]], list[tuple[int, int, str]]]:
    nodes: list[tuple[int, str]] = []
    typed_edges: list[tuple[int, int, str]] = []

    def vertex() -> int:
        result = len(nodes)
        nodes.append((result, "V"))
        return result

    def visit(value: Tree) -> int:
        if not value:
            return vertex()
        left, right = value
        root = vertex()
        left_root = visit(left)
        right_root = visit(right)
        typed_edges.extend([(left_root, root, "L"), (right_root, root, "R")])
        return root

    tree_root = visit(tree)
    exterior_root = vertex()
    typed_edges.append((tree_root, exterior_root, "I"))
    return nodes, typed_edges


def tree_gml(tree: Tree, *, name: str = "tree") -> str:
    """Encode a valid tree as a directed typed incidence graph in GML."""

    nodes, typed_edges = _typed_tree_graph(tree)
    incidence_edges = []
    edge_id = len(nodes)
    for source, target, edge_type in typed_edges:
        nodes.append((edge_id, edge_type))
        incidence_edges.extend(
            [(source, edge_id, "src"), (edge_id, target, "tgt")]
        )
        edge_id += 1
    return f'''graph [
  graphID "{name}"
{_gml_body(nodes, incidence_edges)}
]
'''


def correspondence_certificate() -> dict:
    return {
        "status": "VERIFIED_BY_CONSTRUCTION_AND_DIRECT_TEST",
        "typed_graph": (
            "each directed I/L/R edge is an edge-node with distinct src/tgt "
            "incidence labels; this is fully faithful for finite typed multigraphs"
        ),
        "rule_effect": (
            "delete the selected typed edge, preserve both endpoints, add one "
            "binary vertex, one leaf, and ordered L/R plus replacement T edges"
        ),
        "dpo_sqpo_direct_correspondence": (
            "the only deleted graph item is the selected edge-node and its two "
            "explicitly matched incidence bonds. On valid trees it has no other "
            "incidences, so the DPO dangling condition always holds; SqPO and DPO "
            "therefore produce the same target for every injective match"
        ),
        "jump_closure": (
            "DPO jump closure is I <- K -> I. On valid trees complete deleted "
            "edge-nodes make dangling automatic, also for identity-query "
            "composites. The checked full-input identity reduction has exactly "
            "the same admissible-match count on this language"
        ),
        "guards": (
            "no explicit rule NoEdge is needed on the declared valid-tree language. "
            "The negative PRBF fragment is enforced on every overlap pushout and "
            "on composite input/output; positive totality is discharged by the "
            "recursive valid-tree semantic quotient"
        ),
        "overlap_negative_constraints": TYPED_INCIDENCE_SPECIFICATION,
    }
