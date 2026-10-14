"""Parse the supported MolFoundry jump-observable fragment without model imports."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
from mfobs.core.rules import is_identity
Edge=tuple[int,int,str]
Pair=tuple[int,int]
def _pair(a: int, b: int) -> Pair:
    return (a, b) if a < b else (b, a)

@dataclass(frozen=True)
class ParsedObservable:
    labels: Mapping[int, str]
    edges: tuple[Edge, ...]
    no_edges: tuple[Pair, ...]
    generic_conditions: tuple[str, ...] = ()
    jump_identity: bool = True

    @property
    def edge_set(self) -> set[Edge]:
        return set(self.edges)

def parse_observable(rule) -> ParsedObservable:
    """Extract the input pattern of a MolFoundry DPO jump observable."""

    no_edges = tuple(
        sorted({_pair(c.source, c.target) for c in rule.matchConstraintsNoEdge})
    )
    generic: list[str] = []
    if rule.matchCondition is not None:
        generic.append("matchCondition")
    for key in (
        "matchConstraintsAdjacency",
        "matchConstraintsLabelAny",
        "matchConstraintsShortestPath",
    ):
        if getattr(rule, key):
            generic.append(key)
    return ParsedObservable(
        labels={v.id: v.leftLabel for v in rule.vertices if v.leftLabel is not None},
        edges=tuple(sorted((_pair(e.source.id, e.target.id) + (e.leftLabel,))
                           for e in rule.edges if e.leftLabel is not None)),
        no_edges=no_edges,
        generic_conditions=tuple(generic),
        jump_identity=is_identity(rule),
    )
