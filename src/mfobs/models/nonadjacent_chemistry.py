"""Atom-balanced [2+2] diagnostic with a relational context NoEdge guard.

This is a deliberately selected ISOLATED-ALKENE rate channel, not an
experimentally validated reaction network or an assertion that conjugated
dienes cannot photoreact. In fact, removing only the context guards admits a
valid bicyclobutane product from 1,3-butadiene. The intended diagnostic asks
whether conditions are genuinely active and how pair consumption affects
exact observable closure; no closure or chemical exclusivity is assumed.

All graph labels are ordinary C/H atoms and single/double covalent bonds.
No atom, proton, charge, coordination edge or reporter token is created.
Stereochemistry, wavelength dependence and encounter geometry are omitted.
"""
from __future__ import annotations
from mfobs.core.graph import edge_rows, make_graph, node_labels
from mfobs.core.rules import make_rule
import molfoundry as mf

from fractions import Fraction
from functools import lru_cache


from mfobs.core.gml import disjoint_union
from mfobs.core.molecules import complete_with_h


PRODUCT_NO_EDGES = ((0, 3), (1, 2))
CONTEXT_NO_EDGES = ((0, 2), (1, 3))
PER_MATCH = Fraction(1, 8)


def reaction_sides() -> tuple[mf.Graph, mf.Graph]:
    labels = {u: 'C' for u in range(4)}
    left = make_graph(labels, [(0, 1, '='), (2, 3, '=')])
    right = make_graph(labels, [(0, 1, '-'), (2, 3, '-'), (0, 3, '-'), (1, 2, '-')])
    return left, right


@lru_cache(None)
def guarded_rule() -> mf.Rule:
    """Exclude directly linked alkene pairs as a declared channel choice."""
    return make_rule('isolated-alkene-pair-net-photocycloaddition',
                           *reaction_sides(), PRODUCT_NO_EDGES + CONTEXT_NO_EDGES)


@lru_cache(None)
def ablated_rule() -> mf.Rule:
    """Remove only context guards; never permit parallel product bonds."""
    return make_rule('context-guard-ablated-net-photocycloaddition',
                           *reaction_sides(), PRODUCT_NO_EDGES)


@lru_cache(None)
def ethylene() -> mf.Graph:
    return complete_with_h(make_graph({0: 'C', 1: 'C'}, [(0, 1, '=')]))


@lru_cache(None)
def cyclobutane() -> mf.Graph:
    return complete_with_h(make_graph({u: 'C' for u in range(4)},
                                     [(u, (u+1) % 4, '-') for u in range(4)]))


@lru_cache(None)
def conjugated_diene() -> mf.Graph:
    """1,3-butadiene: 1=0-2=3; context edge 0-2 is actually present."""
    left, _ = reaction_sides()
    return complete_with_h(make_graph(node_labels(left), [*edge_rows(left), (0, 2, '-')]))


@lru_cache(None)
def separated_diene() -> mf.Graph:
    """1,5-hexadiene: 1=0-4-5-2=3, with two saturated spacer carbons."""
    return complete_with_h(make_graph({u: 'C' for u in range(6)},
                                     [(0, 1, '='), (2, 3, '='), (0, 4, '-'),
                                      (4, 5, '-'), (5, 2, '-')]))


@lru_cache(None)
def ethylene_mixture(n: int = 2) -> mf.Graph:
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError('n must be a nonnegative integer')
    return disjoint_union([ethylene()] * n)


def rules() -> dict[str, mf.Rule]:
    return {'isolated_pair': guarded_rule(), 'context_ablated': ablated_rule()}
