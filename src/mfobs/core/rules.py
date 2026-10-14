"""Read native MolFoundry rule spans and construct rules from graph sides."""
from __future__ import annotations

import json
import molfoundry as mf

from mfobs.core.graph import edge_rows, make_graph, node_labels, same_graph


def rule_graph(rule: mf.Rule, side: str = 'left') -> mf.Graph:
    """A native graph for L, K or R, preserving the rule's vertex IDs."""
    if side not in ('left', 'right', 'context'):
        raise ValueError(f'unknown rule side {side!r}')
    nodes, edges = {}, []
    for v in rule.vertices:
        label = v.leftLabel if side != 'right' else v.rightLabel
        if label is not None and (side != 'context' or v.membership.value == 'context'):
            nodes[v.id] = label
    for e in rule.edges:
        label = e.leftLabel if side != 'right' else e.rightLabel
        if label is not None and (side != 'context' or e.isContext):
            edges.append((e.source.id, e.target.id, label))
    return make_graph(nodes, edges)


def rule_vertex_ids(rule: mf.Rule, membership: str) -> tuple[int, ...]:
    return tuple(sorted(v.id for v in rule.vertices if v.membership.value == membership))


def preserves_vertices(rule: mf.Rule) -> bool:
    return all(v.membership.value == 'context' for v in rule.vertices)


def is_identity(rule: mf.Rule) -> bool:
    # A vertex-preserving diagonal span has equal labelled L and R, even
    # when an edge is deleted and recreated across its interface.
    return preserves_vertices(rule) and all(e.leftLabel == e.rightLabel for e in rule.edges)


def removed_edges(rule: mf.Rule) -> set[tuple[int, int, str]]:
    return {(min(e.source.id, e.target.id), max(e.source.id, e.target.id), e.leftLabel)
            for e in rule.edges if e.leftLabel is not None and not e.isContext}


def added_edges(rule: mf.Rule) -> set[tuple[int, int, str]]:
    return {(min(e.source.id, e.target.id), max(e.source.id, e.target.id), e.rightLabel)
            for e in rule.edges if e.rightLabel is not None and not e.isContext}


def restores_input(rule: mf.Rule) -> bool:
    """Whether R restores L, allowing fresh IDs only for recreated vertices."""
    deleted, created = rule_vertex_ids(rule, 'deleted'), rule_vertex_ids(rule, 'created')
    if rule_vertex_ids(rule, 'changed') or len(deleted) != len(created):
        return False
    left, right = rule_graph(rule, 'left'), rule_graph(rule, 'right')
    if not deleted:
        return same_graph(left, right)
    # Pin surviving IDs while comparing the recreated labelled subgraph.
    def pinned(graph, renamed):
        labels = {u: json.dumps(('new', label) if u in renamed else ('kept', u, label))
                  for u, label in node_labels(graph).items()}
        return make_graph(labels, edge_rows(graph))
    return pinned(left, set(deleted)).isomorphism(pinned(right, set(created))) != 0


def make_rule(name: str, left: mf.Graph, right: mf.Graph, no_edges=(), adjacency=()) -> mf.Rule:
    """Build a native DPO rule; equal IDs and labels determine the interface."""
    lnodes, rnodes = node_labels(left), node_labels(right)
    kept = {u for u in lnodes.keys() & rnodes.keys() if lnodes[u] == rnodes[u]}
    common_edges = set(edge_rows(left)) & set(edge_rows(right))
    def body(nodes, edges):
        lines = [f'node [ id {u} label {json.dumps(label)} ]' for u, label in sorted(nodes.items())]
        lines += [f'edge [ source {u} target {v} label {json.dumps(label)} ]'
                  for u, v, label in sorted(edges)]
        return '\n'.join(lines)
    sections = [('left', {u: label for u, label in lnodes.items() if u not in kept}, set(edge_rows(left)) - common_edges),
                ('context', {u: lnodes[u] for u in kept}, common_edges),
                ('right', {u: label for u, label in rnodes.items() if u not in kept}, set(edge_rows(right)) - common_edges)]
    text = 'rule [\nruleID ' + json.dumps(name) + '\n'
    text += ''.join(f'{side} [\n{body(nodes, edges)}\n]\n' for side, nodes, edges in sections)
    text += ''.join(f'constrainNoEdge [ source {u} target {v} ]\n' for u, v in no_edges)
    text += ''.join(c.gml() for c in adjacency)
    return mf.Rule.fromGMLString(text + ']\n', add=False)
