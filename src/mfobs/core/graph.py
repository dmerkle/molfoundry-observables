"""Native MolFoundry graph helpers and concrete applications.

External GML IDs identify vertices across applications. ``graph_key`` compares
ID-labelled states; MolFoundry's ``canonKey`` instead identifies isomorphism
classes. The engine owns graph storage, matching and rule application.
"""
from __future__ import annotations

from collections import Counter
import json
from typing import Iterable, Mapping

import molfoundry as mf

Edge = tuple[int, int, str]


def vertex_ids(graph: mf.Graph) -> dict[int, int]:
    """Internal vertex IDs mapped to unique external IDs, when present."""
    vertices = list(graph.vertices)
    external = [getattr(v, 'externalId', None) for v in vertices]
    use_external = None not in external and len(set(external)) == len(vertices)
    return {v.id: external[i] if use_external else v.id for i, v in enumerate(vertices)}


def node_labels(graph: mf.Graph) -> dict[int, str]:
    ids = vertex_ids(graph)
    return {ids[v.id]: v.stringLabel for v in graph.vertices}


def node_rows(graph: mf.Graph) -> tuple[tuple[int, str], ...]:
    return tuple(sorted(node_labels(graph).items()))


def edge_rows(graph: mf.Graph) -> tuple[Edge, ...]:
    ids = vertex_ids(graph)
    return tuple(sorted((min(ids[e.source.id], ids[e.target.id]),
                         max(ids[e.source.id], ids[e.target.id]), e.stringLabel)
                        for e in graph.edges))


def graph_key(graph: mf.Graph) -> tuple:
    """Exact ID-labelled state key, independent of Python object identity."""
    return node_rows(graph), edge_rows(graph)


def same_graph(left: mf.Graph, right: mf.Graph) -> bool:
    return graph_key(left) == graph_key(right)


def graph_data(graph: mf.Graph) -> dict:
    return {"nodes": [list(row) for row in node_rows(graph)],
            "edges": [list(row) for row in edge_rows(graph)]}


def has_edge(graph: mf.Graph, u: int, v: int, label: str | None = None) -> bool:
    endpoints = min(u, v), max(u, v)
    return any((a, b) == endpoints and (label is None or bond == label)
               for a, b, bond in edge_rows(graph))


def make_graph(nodes: Mapping[int, str] | Iterable[tuple[int, str]],
               edges: Iterable[Edge] = ()) -> mf.Graph:
    """Construct a native graph from labelled vertex/edge rows."""
    nodes = dict(nodes)
    rows = sorted({(min(u, v), max(u, v), label) for u, v, label in edges})
    if any(u == v for u, v, _ in rows):
        raise ValueError("self-loops are outside the simple graph fragment")
    if any(u not in nodes or v not in nodes for u, v, _ in rows):
        raise ValueError("edge endpoint is not a graph node")
    if len({(u, v) for u, v, _ in rows}) != len(rows):
        raise ValueError("parallel edges are outside the simple graph fragment")
    body = [f'node [ id {u} label {json.dumps(label)} ]' for u, label in sorted(nodes.items())]
    body += [f'edge [ source {u} target {v} label {json.dumps(label)} ]' for u, v, label in rows]
    return mf.Graph.fromGMLString('graph [\n' + '\n'.join(body) + '\n]', add=False)


Match = dict[int, int]


def matches(pattern: mf.Graph, host: mf.Graph, no_edges=()) -> list[Match]:
    """MolFoundry injective matches, retaining every map and automorphism."""
    from mfobs.core.molfoundry_execution import matches as engine_matches
    return engine_matches(pattern, host, no_edges)


def apply_match(rule: mf.Rule, host: mf.Graph, match: Match, semantics="DPO") -> mf.Graph | None:
    """Delegate one DPO application, including all conditions, to MolFoundry."""
    if semantics != "DPO":
        raise ValueError("Only MolFoundry DPO application is supported")
    from mfobs.core.molfoundry_execution import apply_match as engine_apply
    return engine_apply(rule, host, match)


def transitions(rule: mf.Rule, host: mf.Graph, semantics="DPO") -> Counter[mf.Graph]:
    """MolFoundry target states with raw admissible-match multiplicities."""
    if semantics != "DPO":
        raise ValueError("Only MolFoundry DPO application is supported")
    from mfobs.core.molfoundry_execution import transitions as engine_transitions
    return engine_transitions(rule, host)


def observable(pattern: mf.Graph, state: mf.Graph, no_edges=()) -> int:
    """Count injective embeddings, the stochastic rewriting observable."""

    from mfobs.core.molfoundry_execution import count_matches
    return count_matches(pattern, state, no_edges)
