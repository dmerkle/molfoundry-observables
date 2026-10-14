"""Finite unary adjacency data and symbolic reductions; execution by MolFoundry.

This module supports conjunctions of local neighbour-count predicates and
NoEdge conditions.  No wildcard/term semantics, nested conditions, shortest-path
conditions, or automatic invariant discovery are implied.  Unknown constructs
raise UnsupportedCondition; the closure reducer must return UNKNOWN for them.
All neighbours are counted in the COMPLETE host, not merely in the matched
positive pattern.  Simple graphs mean neighbours are injectively distinct.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from typing import Mapping
import molfoundry as mf
from mfobs.core.graph import node_labels, edge_rows

class UnsupportedCondition(ValueError):
    pass

OPS = {'<': lambda x,y:x<y, '<=':lambda x,y:x<=y,
       '=':lambda x,y:x==y,
       '>':lambda x,y:x>y, '>=':lambda x,y:x>=y}

@dataclass(frozen=True, order=True)
class Adjacency:
    center: int
    node_labels: tuple[str,...]
    edge_labels: tuple[str,...]
    op: str
    count: int

    def __post_init__(self):
        if self.op not in OPS: raise UnsupportedCondition(f'Unsupported adjacency operator {self.op!r}')
        if not isinstance(self.count,int) or isinstance(self.count,bool):
            raise UnsupportedCondition('Adjacency count must be an integer')
        # Symbolic counting requires explicit nonempty finite label filters.
        if not self.node_labels or not self.edge_labels:
            raise UnsupportedCondition('Explicit conditions require explicit nonempty node/edge label filters')
        if any(not isinstance(s,str) for s in (*self.node_labels,*self.edge_labels)):
            raise UnsupportedCondition('Only string label filters are supported')

    def reindex(self, mapping: Mapping[int,int]) -> 'Adjacency':
        return Adjacency(mapping[self.center],self.node_labels,self.edge_labels,self.op,self.count)

    def accepts_count(self, n:int) -> bool:
        return bool(OPS[self.op](n,self.count))

    def count_at(self, host:mf.Graph, image:int) -> int:
        labels=node_labels(host)
        return sum(bond in self.edge_labels and
                   labels[v if u==image else u] in self.node_labels
                   for u,v,bond in edge_rows(host) if image in (u,v))

    def accepts(self, host:mf.Graph, image:int) -> bool:
        return self.accepts_count(self.count_at(host,image))

    def gml(self) -> str:
        quote=lambda s:json.dumps(s,ensure_ascii=True)
        ns=' '.join('label '+quote(s) for s in self.node_labels)
        es=' '.join('label '+quote(s) for s in self.edge_labels)
        return (f'  constrainAdj [\n    id {self.center}\n    op {quote(self.op)}\n'
                f'    count {self.count}\n    nodeLabels [ {ns} ]\n'
                f'    edgeLabels [ {es} ]\n  ]\n')

    @classmethod
    def from_constraint(cls, constraint) -> 'Adjacency':
        """From MolFoundry's ``AdjacencyConstraint``; an omitted label filter fails closed."""
        labels=lambda value:tuple(sorted(set(value or ())))
        return cls(constraint.id,labels(constraint.nodeLabels),labels(constraint.edgeLabels),
                   constraint.op.value,constraint.count)

    def engine_constraint(self):
        """As MolFoundry's ``AdjacencyConstraint``."""
        import molfoundry as mf
        return mf.AdjacencyConstraint(self.center, mf.ConstraintOp(self.op), self.count,
                                      nodeLabels=self.node_labels, edgeLabels=self.edge_labels)

def identity_rule(graph, name='identity', no_edges=(), adjacency=()):
    """MolFoundry's identity rule ``G <- G -> G`` on ``graph`` (``Rule.identity``)
    under NoEdge pairs and adjacency conditions: the observable counting the
    pattern's admissible matches. Rule vertices keep the graph's node IDs."""
    import molfoundry as mf
    constraints=[mf.NoEdgeConstraint(u,v) for u,v in no_edges]
    constraints+=[c.engine_constraint() for c in adjacency]
    return mf.Rule.identity(graph, name,
                            constraints=constraints, add=False)


def matches(pattern, host, no_edges=(), adjacency=()):
    from mfobs.core.molfoundry_execution import matches as engine_matches
    return engine_matches(pattern, host, no_edges, adjacency)


def transitions(rule, host):
    from mfobs.core.molfoundry_execution import transitions as engine_transitions
    return engine_transitions(rule, host)


def engine_has_conditions(rule):
    """Whether a MolFoundry rule carries any match constraint or application condition."""
    return rule.matchCondition is not None or any(getattr(rule,key) for key in (
        'matchConstraintsNoEdge','matchConstraintsAdjacency',
        'matchConstraintsLabelAny','matchConstraintsShortestPath'))


def engine_no_edges(rule):
    """NoEdge pairs of a MolFoundry rule, including those implied by created edges."""
    return tuple((c.source,c.target) for c in rule.matchConstraintsNoEdge)


def engine_conditions(rule):
    """Read only the supported condition fragment; never silently drop others."""
    if rule.matchCondition is not None:
        raise UnsupportedCondition('Nested/generic matchCondition is not supported')
    for key in ['matchConstraintsLabelAny','matchConstraintsShortestPath']:
        if getattr(rule,key):raise UnsupportedCondition(f'{key} is not supported')
    return tuple(Adjacency.from_constraint(c) for c in rule.matchConstraintsAdjacency)


def adjacency_contradictions(rule):
    """Return certified contradictions in the supported adjacency fragment.

    Conditions with the same centre and filters constrain one nonnegative
    integer.  Positive pattern edges also provide a lower bound.  This is a
    deliberately small exact solver, not a general nested-condition prover.
    """

    from mfobs.core.observables import parse_observable as base_parse

    parsed = base_parse(rule)
    conditions = engine_conditions(rule)
    grouped = {}
    for condition in conditions:
        key = (condition.center, condition.node_labels, condition.edge_labels)
        grouped.setdefault(key, []).append(condition)
    contradictions = []
    labels = parsed.labels
    for (center, neighbour_labels, edge_labels), values in grouped.items():
        lower = 0
        upper = None
        for condition in values:
            if condition.op == '<':
                candidate = condition.count - 1
                upper = candidate if upper is None else min(upper, candidate)
            elif condition.op == '<=':
                upper = condition.count if upper is None else min(upper, condition.count)
            elif condition.op == '=':
                lower = max(lower, condition.count)
                upper = condition.count if upper is None else min(upper, condition.count)
            elif condition.op == '>':
                lower = max(lower, condition.count + 1)
            elif condition.op == '>=':
                lower = max(lower, condition.count)
        positive = sum(
            bond in edge_labels
            and labels[other] in neighbour_labels
            for left, right, bond in parsed.edges
            if center in (left, right)
            for other in (right if left == center else left,)
        )
        lower = max(lower, positive)
        if upper is not None and lower > upper:
            contradictions.append({
                'center': center,
                'node_labels': list(neighbour_labels),
                'edge_labels': list(edge_labels),
                'lower_bound': lower,
                'upper_bound': upper,
                'reason': 'empty-nonnegative-integer-interval',
            })
    return tuple(contradictions)


def parse_observable(rule):
    from dataclasses import replace
    from mfobs.core.observables import parse_observable as base_parse
    p=base_parse(rule)
    adj=engine_conditions(rule)
    if any(c.center not in p.labels for c in adj):
        raise UnsupportedCondition('Transported adjacency center outside observable input')
    p=replace(p,generic_conditions=tuple(c for c in p.generic_conditions if c!='matchConstraintsAdjacency'))
    return p,adj
