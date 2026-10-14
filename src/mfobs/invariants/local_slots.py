"""Symbolic local-slot reduction, with NO concrete states or orbit evaluations.

This is an explicitly supported invariant fragment, not an invariant discoverer
or a complete solver for arbitrary graph constraints. Exact results preserve raw
injective embedding counts. Unsupported partial ownership stays UNKNOWN.
"""
from __future__ import annotations
from mfobs.core.graph import node_rows, make_graph, edge_rows, has_edge, node_labels
from mfobs.core.rules import make_rule, rule_graph, preserves_vertices
from collections import Counter, defaultdict
from dataclasses import dataclass
import json

from mfobs.core.conditions import UnsupportedCondition, engine_no_edges, parse_observable


CONTRACT = {
    'kind': 'local-owned-slots-v1',
    'node_labels': ['C', 'P', 'L'],
    'edge_types': [['C', 'P', 'cp'], ['C', 'L', 'cl'], ['C', 'C', 'next'], ['P', 'L', 'bound']],
    'each_C_exactly_one_P_and_one_L': True,
    'each_P_L_at_most_one_C': True,
    'each_P_L_at_most_one_bound_partner': True,
    'bound_endpoints_have_the_same_C': True,
    'next_max_degree': None,
}


@dataclass
class Reduction:
    status: str
    reason: str
    key: str | None = None
    rule: object = None
    quotient: dict | None = None
    removed_no_edges: tuple = ()

    def record(self):
        return {'status': self.status, 'reason': self.reason, 'key': self.key,
                'gml': self.rule.getGMLString() if self.rule is not None else None,
                'quotient': self.quotient, 'removed_no_edges': list(self.removed_no_edges),
                'multiplicity_factor': 1 if self.status == 'exact' else None}


class Normalizer:
    def __init__(self, contract=None):
        self.contract = dict(CONTRACT if contract is None else contract)
        base = {k: v for k, v in self.contract.items() if k != 'next_max_degree'}
        if base != {k: v for k, v in CONTRACT.items() if k != 'next_max_degree'}:
            raise ValueError('Unsupported local-slot invariant fragment')
        degree = self.contract.get('next_max_degree')
        if degree is not None and (type(degree) is not int or degree < 0):
            raise ValueError('next_max_degree must be a nonnegative integer or null')
        self.next_max_degree = degree
        self.cache = {}
        self.stats = Counter()

    def _structure(self, graph):
        """Positive contradictions are zero. Missing ownership is not zero."""
        labels = node_labels(graph)
        if set(labels.values()) - {'C', 'P', 'L'}:
            return 'zero', 'node label forbidden by local signature', None
        by_center = {u: {'P': set(), 'L': set()} for u, l in node_rows(graph) if l == 'C'}
        owners = {u: set() for u, l in node_rows(graph) if l in ('P', 'L')}
        bound = defaultdict(set); following = defaultdict(set)
        for u, v, edge in edge_rows(graph):
            a, b = labels[u], labels[v]
            if edge in ('cp', 'cl'):
                partner = 'P' if edge == 'cp' else 'L'
                if {a, b} != {'C', partner}:
                    return 'zero', 'ownership edge has forbidden endpoint types', None
                c, p = (u, v) if a == 'C' else (v, u)
                by_center[c][partner].add(p); owners[p].add(c)
            elif edge == 'next' and a == b == 'C':
                following[u].add(v); following[v].add(u)
            elif edge == 'bound' and {a, b} == {'P', 'L'}:
                bound[u].add(v); bound[v].add(u)
            else:
                return 'zero', 'edge label/endpoint types forbidden by local signature', None
        if any(len(v) > 1 for slot in by_center.values() for v in slot.values()):
            return 'zero', 'one C has two distinct P or L vertices', None
        if any(len(v) > 1 for v in owners.values()):
            return 'zero', 'one P/L has two distinct C owners', None
        if any(len(v) > 1 for v in bound.values()):
            return 'zero', 'one P/L has two distinct bound partners', None
        if self.next_max_degree is not None and any(len(v) > self.next_max_degree for v in following.values()):
            return 'zero', 'positive next-degree exceeds the supplied bound', None
        # Known different owners already contradict the same-owner bound axiom.
        for u, v, edge in edge_rows(graph):
            if edge == 'bound' and owners[u] and owners[v] and owners[u] != owners[v]:
                return 'zero', 'bound edge joins different known owners', None
        if any(not v for v in owners.values()) or any(not v for slot in by_center.values() for v in slot.values()):
            return 'unknown', 'partial ownership: this normalizer does not guess or merge missing slot vertices', None
        return 'exact', 'complete anchored slots satisfy all positive local bounds', (by_center, owners, following)

    def normalize_graph(self, graph, no_edges=(), adjacency=()):
        if adjacency:
            return Reduction('unknown', 'adjacency predicates are outside this first local-slot normalizer')
        labels = node_labels(graph)
        if any(len(pair) != 2 or any(u not in labels for u in pair) for pair in no_edges):
            return Reduction('unknown', 'invalid NoEdge endpoints')
        if any(a == b for a, b in no_edges):
            return Reduction('unknown', 'self NoEdge not supported')
        if any(has_edge(graph, u, v) for u, v in no_edges):
            return Reduction('zero', 'positive edge contradicts NoEdge')
        status, reason, structure = self._structure(graph)
        if status != 'exact':
            return Reduction(status, reason)
        slots, owners, following = structure
        retained = set(); removed = []
        for u, v in no_edges:
            a, b = labels[u], labels[v]
            if a == b == 'C':
                retained.add(tuple(sorted((u, v))))
            elif {a, b} == {'P', 'L'} and owners[u] == owners[v]:
                retained.add(tuple(sorted((u, v))))
            else:
                # For complete anchored triples the only other possible typed
                # edges are already positive ownership edges (checked above).
                removed.append((u, v))
        centers = sorted(slots); index = {u: i for i, u in enumerate(centers)}
        states = []
        for c in centers:
            p, l = next(iter(slots[c]['P'])), next(iter(slots[c]['L']))
            states.append(1 if has_edge(graph, p, l, 'bound') else
                          0 if tuple(sorted((p, l))) in retained else None)
        quotient = {'occupancies': states,
                    'next_edges': sorted([index[u], index[v]] for u, v, label in edge_rows(graph) if label == 'next'),
                    'no_next_edges': sorted([index[u], index[v]] for u, v in retained if labels[u] == labels[v] == 'C'),
                    'source_center_ids': centers, 'site_count': len(centers)}
        identity = make_rule('locally-normalized-count', graph, graph, tuple(sorted(retained)))
        rule = identity
        return Reduction('exact', 'same input pattern; only locally implied NoEdges removed; vertex-preserving jump flattened',
                         rule.canonKey, rule, quotient, tuple(sorted(removed)))

    def normalize(self, rule):
        key = rule.canonKey
        if key not in self.cache:
            try:
                p, adjacency = parse_observable(rule)
                if not p.jump_identity or p.generic_conditions or not preserves_vertices(rule):
                    value = Reduction('unknown', 'only vertex-preserving supported diagonal spans are certified')
                else:
                    value = self.normalize_graph(make_graph(p.labels, p.edges), p.no_edges, adjacency)
            except (UnsupportedCondition, ValueError, KeyError) as exc:
                value = Reduction('unknown', str(exc))
            self.cache[key] = value
            self.stats[value.status] += 1
        return self.cache[key]

    def certify_rule_preservation(self, rule):
        """Check sufficient local preservation premises; no host enumeration."""
        left, right = rule_graph(rule), rule_graph(rule, "right")
        if node_labels(left) != node_labels(right) or not preserves_vertices(rule):
            raise ValueError('Preservation checker requires all vertices to be preserved')
        if {e for e in edge_rows(left) if e[2] != 'bound'} != {e for e in edge_rows(right) if e[2] != 'bound'}:
            raise ValueError('Ownership/next topology changes are not covered')
        for graph in (left, right):
            status, reason, _ = self._structure(graph)
            if status != 'exact':
                raise ValueError('Unsupported rule-side structure: '+reason)
        no_edges = {tuple(sorted(pair)) for pair in engine_no_edges(rule)}
        added, deleted = set(edge_rows(right))-set(edge_rows(left)), set(edge_rows(left))-set(edge_rows(right))
        if any(edge != 'bound' or tuple(sorted((u, v))) not in no_edges for u, v, edge in added):
            raise ValueError('Added bound edge requires an explicit absence guard')
        if any(edge != 'bound' for _, _, edge in deleted):
            raise ValueError('Only bound-edge deletions covered')
        return {'passed': True, 'vertices_preserved': True, 'ownership_and_next_preserved': True,
                'added_bound_edges': list(added), 'deleted_bound_edges': list(deleted),
                'argument': 'Unique owned P/L endpoints have no other permitted bound partner. NoEdge permits adding their sole optional edge; deletion preserves all upper bounds. Read-only cooperative context does not change topology.',
                'method': 'syntactic sufficient-condition check, not sampled-state testing'}


def meaning(quotient):
    values = quotient['occupancies']; n = len(values)
    if not n: return 'empty-pattern constant 1'
    if n == 1: return {0: '0 (free)', 1: '1 (occupied)', None: '* (unspecified)'}[values[0]]
    adj = {i: set() for i in range(n)}
    for u, v in quotient['next_edges']: adj[u].add(v); adj[v].add(u)
    ends = [u for u in adj if len(adj[u]) == 1]
    if len(ends) == 2 and all(len(v) <= 2 for v in adj.values()) and not quotient['no_next_edges']:
        order = [ends[0]]
        while len(order) < n:
            following = [v for v in adj[order[-1]] if v not in order]
            if len(following) != 1: break
            order.append(following[0])
        if len(order) == n:
            word = ''.join('*' if values[u] is None else str(values[u]) for u in order)
            return min(word, word[::-1])+' (path motif)'
    return f'{n}-site conditional motif: '+json.dumps({k: quotient[k] for k in ('occupancies','next_edges','no_next_edges')}, separators=(',', ':'))
