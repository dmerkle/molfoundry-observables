"""Aromatic-ring chlorination: Ω, the neighbour-Cl NAC, growth and selected C-Cl count."""
from __future__ import annotations
from mfobs.core.graph import edge_rows, make_graph, node_rows
from mfobs.core.rules import make_rule
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial
from mfobs.core.conditions import Adjacency

from mfobs.core.gml import disjoint_union
from mfobs.state_space import StateSpace, covalent_theory
import molfoundry as mf
from mfobs.core.pattern_support import Pattern


def state_family(chlorine_total=4):
    """Declare Ω for this model, independently of its selected queries."""
    if type(chlorine_total) is not int or chlorine_total < 4:
        raise ValueError("chlorine_total must be an integer >= 4")
    return StateSpace("chlorinated-rings", {
        "ring_vertices": 5, "substitution_sites": 4,
        "head": "one alkyl-substituted cyclopentadienyl anion with aromatic bonds",
        "tail": "one nonempty finite linear alkyl tail with explicit hydrogens",
        "chlorine_total": chlorine_total,
        "HCl": "chlorine_total minus ring-bound Cl",
        "pools": {"methane": "arbitrary finite", "H2": "arbitrary finite"},
    }, covalent_theory())


OMEGA = state_family()


RING_VERTICES = OMEGA.definition["ring_vertices"]


SUBSTITUTION_SITES = OMEGA.definition["substitution_sites"]


BOND_ORDERS = {'-': Fraction(1), ':': Fraction(3, 2)}


# Raw labelled-H matches: 3! on the terminal group × 4! on methane.
RATE_FACTORS = {'grow': Fraction(1, factorial(3) * factorial(4))}


def bond_order_sums(graph: mf.Graph) -> Counter:
    """Conserved graph weights; NOT chemical valence for aromatic C-."""
    result = Counter()
    for u, v, label in edge_rows(graph):
        result[u] += BOND_ORDERS[label]
        result[v] += BOND_ORDERS[label]
    return result


def formula(graph: mf.Graph) -> dict[str, int]:
    """Element counts, with formal charge reported separately."""
    return dict(sorted(Counter('C' if label == 'C-' else label
                               for _, label in node_rows(graph)).items()))


def formal_charge(graph: mf.Graph) -> int:
    return -sum(label == 'C-' for _, label in node_rows(graph))


def _validate(n, mask):
    if type(n) is not int or n < 0:
        raise ValueError('n must be a nonnegative integer')
    if type(mask) is not int or not 0 <= mask < 1 << SUBSTITUTION_SITES:
        raise ValueError('mask must be an integer in 0..15')


@lru_cache(None, typed=True)
def molecule(n: int = 0, mask: int = 0) -> mf.Graph:
    """C(6+n) H(7+2n-popcount(mask)) Cl(popcount(mask)), charge -1.

    Ring atoms 0..4, charged anchor 0, tail atoms 5..5+n. Mask bits
    correspond to neutral ring sites 1..4. Hydrogen completion is explicit
    and structure-aware; no neutral-carbon valence heuristic is used for C-.
    """
    _validate(n, mask)
    nodes = {u: 'C' for u in range(6+n)}
    nodes[0] = 'C-'
    edges = [(u, (u+1) % RING_VERTICES, ':') for u in range(RING_VERTICES)]
    edges += [(0, 5, '-')] + [(u, u+1, '-') for u in range(5, 5+n)]
    nxt = 6+n
    for position in range(1, SUBSTITUTION_SITES + 1):
        nodes[nxt] = 'Cl' if mask & (1 << (position-1)) else 'H'
        edges.append((position, nxt, '-'))
        nxt += 1
    for position in range(5, 6+n):
        for _ in range(3 if position == 5+n else 2):
            nodes[nxt] = 'H'
            edges.append((position, nxt, '-'))
            nxt += 1
    result = make_graph(nodes, edges)
    OMEGA.validate_graph(result)
    return result


@lru_cache(None)
def methane() -> mf.Graph:
    return make_graph({0: 'C', 1: 'H', 2: 'H', 3: 'H', 4: 'H'},
                      [(0, i, '-') for i in range(1, 5)])


@lru_cache(None)
def hcl() -> mf.Graph:
    return make_graph({0: 'H', 1: 'Cl'}, [(0, 1, '-')])


@lru_cache(None)
def h2() -> mf.Graph:
    return make_graph({0: 'H', 1: 'H'}, [(0, 1, '-')])


_REAGENT_BUILDERS = {'methane': methane, 'hcl': hcl, 'h2': h2}


def host(n=0, mask=0, methane=0, hcl=0, h2=0) -> mf.Graph:
    for amount in (methane, hcl, h2):
        if type(amount) is not int or amount < 0:
            raise ValueError('Reagent/product populations must be nonnegative integers')
    return disjoint_union([molecule(n, mask)] + [_REAGENT_BUILDERS['methane']()]*methane
                          + [_REAGENT_BUILDERS['hcl']()]*hcl
                          + [_REAGENT_BUILDERS['h2']()]*h2)


@lru_cache(None)
def grow_rule() -> mf.Rule:
    nodes = {0: 'C', 1: 'H', 2: 'C', 3: 'H', 4: 'H', 5: 'H',
             6: 'H', 7: 'H', 8: 'H'}
    left_edges = [(0, 1, '-'), (0, 7, '-'), (0, 8, '-'),
                  (2, 3, '-'), (2, 4, '-'), (2, 5, '-'), (2, 6, '-')]
    right_edges = [edge for edge in left_edges
                   if edge not in ((0, 1, '-'), (2, 3, '-'))]
    right_edges += [(0, 2, '-'), (1, 3, '-')]
    # The already attached neighbour may be charged Cp C- at n=0 or a normal
    # tail C later. Omitting a positive neighbour vertex avoids fake atom types.
    attached_terminal = Adjacency(0, ('C', 'C-'), ('-',), '=', 1)
    return make_rule('cp-terminal-methane-homologation-net-channel',
                           make_graph(nodes, left_edges), make_graph(nodes, right_edges),
                           (), (attached_terminal,))


def seed_graph() -> mf.Graph:
    return make_graph({0: 'C', 1: 'Cl'}, [(0, 1, '-')])


RULE_BODY = r"""
rule [
  ruleID "single-chlorination-neighbour-Cl-NAC"
  left [
    edge [ source 0 target 1 label "-" ]
    edge [ source 2 target 3 label "-" ]
  ]
  context [
    node [ id 0 label "C" ]
    node [ id 1 label "H" ]
    node [ id 2 label "H" ]
    node [ id 3 label "Cl" ]
  ]
  right [
    edge [ source 0 target 3 label "-" ]
    edge [ source 1 target 2 label "-" ]
  ]
  constrainAdj [
    id 0
    op "="
    count 2
    nodeLabels [ label "C" label "C-" ]
    edgeLabels [ label ":" ]
  ]
"""

NEIGHBOUR_NAC = r"""
  condition [
    not [
      exists [
        node [ id 4 label "C" ]
        node [ id 5 label "Cl" ]
        edge [ source 0 target 4 label ":" ]
        edge [ source 4 target 5 label "-" ]
      ]
    ]
  ]
"""

RULE_GML = RULE_BODY + NEIGHBOUR_NAC + "]\n"


def chlorination_rule(guarded=True):
    """The same native rule, with or without the neighbour-Cl NAC."""
    text = RULE_GML if guarded else RULE_BODY + "]\n"
    return mf.Rule.fromGMLString(text, add=False)


def _identity(name: str, body: str):
    graph = mf.Graph.fromGMLString(f"graph [ {body} ]", add=False)
    return mf.Rule.identity(graph, name, add=False)


def query_patterns():
    """The selected query; no auxiliary counting patterns are supplied."""

    return {
        "Z": Pattern(
            "Z",
            _identity(
                "Z",
                'node [ id 0 label "C" ] node [ id 1 label "Cl" ] '
                'edge [ source 0 target 1 label "-" ]',
            ),
            "positive",
            {"meaning": "chlorinated neutral ring site"},
        ),
    }
