"""Matching and DPO application through MolFoundry's public Python API.

Counts retain every injective map. Applications preserve the host's external
vertex IDs; created vertices receive the smallest unused IDs. MolFoundry
validates each supplied match, including edited mappings.
"""
from collections import Counter

import molfoundry as mf

from mfobs.core.graph import edge_rows, graph_key, make_graph, node_labels, vertex_ids
from mfobs.core.rules import rule_vertex_ids


def canonical_key(graph: mf.Graph):
    return str(graph.canonKey)


def _counting_rule(pattern, no_edges, adjacency):
    if isinstance(pattern, mf.Rule):
        if no_edges or adjacency:
            raise ValueError('A native rule already contains its conditions')
        return pattern
    from mfobs.core.conditions import identity_rule
    return identity_rule(pattern, 'count', tuple(no_edges), tuple(adjacency))


def count_matches(pattern, host: mf.Graph, no_edges=(), adjacency=()):
    """Raw admissible matches in this concrete host, without extra host copies."""
    rule = _counting_rule(pattern, no_edges, adjacency)
    return rule.numMatches([host], mf.ApplyMode.UseAny)


class Match(dict):
    """An enumerated match; preserve engine failures for an unchanged mapping."""
    __slots__ = ('enumerated_for',)


def matches(pattern, host: mf.Graph, no_edges=(), adjacency=()):
    """Every rule/pattern match as rule vertex ID -> host vertex ID."""
    rule = _counting_rule(pattern, no_edges, adjacency)
    ids, provenance = vertex_ids(host), (rule.id, graph_key(host))
    result = []
    for row in rule.matches([host], mf.ApplyMode.UseAny):
        match = Match({node: ids[vertex] for node, (_, vertex) in row.items()})
        match.enumerated_for = (*provenance, frozenset(match.items()))
        result.append(match)
    return result


def apply_match(rule: mf.Rule, host: mf.Graph, match):
    """Apply at the supplied match; return the native product, or None if invalid."""
    left_ids = {v.id for v in rule.vertices if v.leftLabel is not None}
    if (set(match) != left_ids or len(set(match.values())) != len(match)
            or not set(match.values()) <= set(node_labels(host))):
        return None
    validated = getattr(match, 'enumerated_for', None) == (
        rule.id, graph_key(host), frozenset(match.items()))
    host_ids = vertex_ids(host)
    internal = {external: vertex for vertex, external in host_ids.items()}
    spec = {node: (0, internal[target]) for node, target in match.items()}
    with mf.DG(graphDatabase=[host]).build() as builder:
        try:
            _, vertex_map = builder.applyMatch([host], rule, spec)
        except mf.LogicError:
            if validated:
                raise
            return None
    created, used, fresh = {}, set(node_labels(host)), 0
    for node in rule_vertex_ids(rule, 'created'):
        while fresh in used:
            fresh += 1
        created[node] = fresh
        used.add(fresh)
    ids, nodes, products = {}, node_labels(host) if not match else {}, {}
    for entry in vertex_map:
        if entry.target is None:
            continue
        key = (entry.target.graph.id, entry.targetInstance, entry.target.id)
        ids[key] = host_ids[entry.source.id] if entry.source is not None else created[entry.ruleNodeId]
        nodes[ids[key]] = entry.target.stringLabel
        products.setdefault(key[:2], entry.target.graph)
    edges = list(edge_rows(host)) if not match else []
    edges += [(ids[(*copy, e.source.id)], ids[(*copy, e.target.id)], e.stringLabel)
              for copy, product in products.items() for e in product.edges]
    return make_graph(nodes, edges)


def transitions(rule: mf.Rule, host: mf.Graph):
    """Native target states grouped by exact IDs, with raw match multiplicities."""
    representatives, result = {}, Counter()
    before = graph_key(host)
    for match in matches(rule, host):
        target = apply_match(rule, host, match)
        if target is not None and (key := graph_key(target)) != before:
            representative = representatives.setdefault(key, target)
            result[representative] += 1
    return result
