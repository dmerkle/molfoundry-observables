"""Artificial but chemical-bond-only structural NAC oscillator.

The state language is an unbounded family of oriented macrocycles.  Every
repeat contains one pendant carbon ring-shift switch in topology Q0, Q1, or
Q2.  Direction is encoded by a covalent C-S-O-C backbone; switch identity is
encoded by one of three ordinary C-C closure bonds.  No atom or edge is a
state label.

Neighbour-assisted rules implement Q1+Q0->Q1+Q1, Q2+Q1->Q2+Q2, and
Q0+Q2->Q0+Q0 along the backbone.  A second copy of every channel inspects the
next repeat through a genuine NoEdge condition.  Six symmetric neighbour-swap
rules provide tunable mixing without changing the one-site balance.  The
guarded rules are the source of damping in the one-site closure.

Three further condensation rules consume a valence-complete capped precursor,
release mapped water, and insert a copy of the local Q state.  At equal
per-site rate this grows N while cancelling from normalized one-site fractions;
precursor activity is represented by an effective chemostatted rate.

The transformations are hypothetical fuel-driven net ring shifts, not proposed
elementary mechanisms or a claim about a synthesised material.
"""
from __future__ import annotations
from mfobs.core.graph import edge_rows, make_graph, node_labels, node_rows
from mfobs.core.rules import make_rule, rule_graph
import molfoundry as mf

from fractions import Fraction
from functools import lru_cache
import math


from mfobs.core.gml import bond_order_sums, union
from mfobs.core.molecules import VALENCE, complete_with_h, formula


BASE_10 = "base_10"
BASE_21 = "base_21"
BASE_02 = "base_02"
GUARD_10 = "guard_10_not1"
GUARD_21 = "guard_21_not2"
GUARD_02 = "guard_02_not0"
SWAP_01 = "swap_01_to_10"
SWAP_10 = "swap_10_to_01"
SWAP_12 = "swap_12_to_21"
SWAP_21 = "swap_21_to_12"
SWAP_20 = "swap_20_to_02"
SWAP_02 = "swap_02_to_20"
GROW_0 = "grow_copy_0"
GROW_1 = "grow_copy_1"
GROW_2 = "grow_copy_2"
CORE_REACTIONS = (BASE_10, BASE_21, BASE_02, GUARD_10, GUARD_21, GUARD_02)
SWAP_REACTIONS = (SWAP_01, SWAP_10, SWAP_12, SWAP_21, SWAP_20, SWAP_02)
GROW_REACTIONS = (GROW_0, GROW_1, GROW_2)
REACTION_ORDER = (*CORE_REACTIONS, *SWAP_REACTIONS, *GROW_REACTIONS)
RATE_NAMES = {
    BASE_10: "r", BASE_21: "r", BASE_02: "r",
    GUARD_10: "g", GUARD_21: "g", GUARD_02: "g",
    **{name: "m" for name in SWAP_REACTIONS},
    **{name: "beta" for name in GROW_REACTIONS},
}
GROW_STATES = {GROW_0: 0, GROW_1: 1, GROW_2: 2}
SWAP_TRANSITIONS = {
    SWAP_01: ((0, 1), (1, 0)),
    SWAP_10: ((1, 0), (0, 1)),
    SWAP_12: ((1, 2), (2, 1)),
    SWAP_21: ((2, 1), (1, 2)),
    SWAP_20: ((2, 0), (0, 2)),
    SWAP_02: ((0, 2), (2, 0)),
}
RULE_TRANSITIONS = {
    BASE_10: (1, 0, 1, None),
    BASE_21: (2, 1, 2, None),
    BASE_02: (0, 2, 0, None),
    GUARD_10: (1, 0, 1, 1),
    GUARD_21: (2, 1, 2, 2),
    GUARD_02: (0, 2, 0, 0),
}

# Local site IDs are p, c0..c4, and the carbonyl oxygen.  The switch edge
# advances along the carbon chain but remains an ordinary single bond.
P = 0
CORE = (1, 2, 3, 4, 5)
OXO = 6
STATE_EDGES = {
    0: (CORE[0], CORE[2], "-"),
    1: (CORE[1], CORE[3], "-"),
    2: (CORE[2], CORE[4], "-"),
}


def _shift_edge(edge, offset):
    return edge[0] + offset, edge[1] + offset, edge[2]


def _site_nodes(offset=0):
    return {
        offset + P: "C",
        **{offset + node: "C" for node in CORE},
        offset + OXO: "O",
    }


def site_heavy(state: int | None, offset=0) -> mf.Graph:
    """One pendant switch scaffold, optionally with its state closure."""
    if state not in (None, 0, 1, 2):
        raise ValueError("switch state must be 0, 1, 2, or None")
    edges = [(offset + P, offset + CORE[0], "-")]
    edges += [
        (offset + left, offset + right, "-")
        for left, right in zip(CORE, CORE[1:])
    ]
    edges.append((offset + CORE[-1], offset + OXO, "="))
    if state is not None:
        edges.append(_shift_edge(STATE_EDGES[state], offset))
    return make_graph(_site_nodes(offset), edges)


def _site_offset(index: int) -> int:
    return 20 * index


def _connector_ids(index: int) -> tuple[int, int]:
    base = 20 * index
    return base + 10, base + 11


@lru_cache(None)
def capped_site(state: int) -> mf.Graph:
    """Return one complete molecular repeat with valence-completing O-H caps.

    This is a depiction molecule, not an extra state of the macrocycle
    language.  It exposes the same local switch topology without presenting a
    partial query graph as though it were a chemically complete molecule.
    """
    site = site_heavy(state)
    nodes = {**node_labels(site), 10: "O", 11: "S", 12: "O"}
    edges = (
        *edge_rows(site),
        (10, P, "-"),
        (P, 11, "-"),
        (11, 12, "-"),
    )
    return complete_with_h(make_graph(nodes, edges))


@lru_cache(None)
def macrocycle(states: tuple[int, ...]) -> mf.Graph:
    """Return a complete explicit-H cyclic polymer for a state word."""
    states = tuple(states)
    if len(states) < 3 or any(state not in (0, 1, 2) for state in states):
        raise ValueError("a macrocycle needs at least three Q0/Q1/Q2 repeats")
    nodes, edges = {}, []
    for index, state in enumerate(states):
        site = site_heavy(state, _site_offset(index))
        nodes.update(node_labels(site))
        edges.extend(edge_rows(site))
        sulfur, oxygen = _connector_ids(index)
        nodes[sulfur], nodes[oxygen] = "S", "O"
        current = _site_offset(index) + P
        following = _site_offset((index + 1) % len(states)) + P
        edges.extend(
            ((current, sulfur, "-"), (sulfur, oxygen, "-"),
             (oxygen, following, "-"))
        )
    return complete_with_h(make_graph(nodes, edges))


def word_pattern(states: tuple[int, ...]) -> mf.Graph:
    """Oriented connected observable for a consecutive switch word."""
    states = tuple(states)
    if not states or any(state not in (0, 1, 2) for state in states):
        raise ValueError("word pattern needs one or more Q0/Q1/Q2 states")
    nodes, edges = {}, []
    for index, state in enumerate(states):
        site = site_heavy(state, _site_offset(index))
        nodes.update(node_labels(site))
        edges.extend(edge_rows(site))
        if index + 1 < len(states):
            sulfur, oxygen = _connector_ids(index)
            nodes[sulfur], nodes[oxygen] = "S", "O"
            edges.extend(
                ((_site_offset(index) + P, sulfur, "-"),
                 (sulfur, oxygen, "-"),
                 (oxygen, _site_offset(index + 1) + P, "-"))
            )
    # Boundary stubs make the path orientation and word origin explicit.
    incoming_o, outgoing_s = 20 * len(states) + 12, 20 * len(states) + 13
    nodes[incoming_o], nodes[outgoing_s] = "O", "S"
    edges.extend(
        ((incoming_o, P, "-"),
         (_site_offset(len(states) - 1) + P, outgoing_s, "-"))
    )
    return make_graph(nodes, edges)


def observable_rule(states: tuple[int, ...]) -> mf.Rule:
    graph = word_pattern(states)
    name = "observable-Q" + "Q".join(str(state) for state in states)
    return make_rule(name, graph, graph)


def _connector(left_p: int, right_p: int, start: int) -> mf.Graph:
    return make_graph(
        {left_p: "C", right_p: "C", start: "S", start + 1: "O"},
        ((left_p, start, "-"), (start, start + 1, "-"),
         (start + 1, right_p, "-")),
    )


def _merge_overlapping(*graphs: mf.Graph) -> mf.Graph:
    nodes, edges = {}, []
    for graph in graphs:
        for node, label in node_rows(graph):
            if node in nodes and nodes[node] != label:
                raise ValueError("incompatible overlapping node labels")
            nodes[node] = label
        edges.extend(edge_rows(graph))
    return make_graph(nodes, edges)


def _transition_hydrogens(old: int, new: int, offset: int, first_h: int):
    old_endpoints = set(STATE_EDGES[old][:2])
    new_endpoints = set(STATE_EDGES[new][:2])
    gains_h = sorted(old_endpoints - new_endpoints)
    loses_h = sorted(new_endpoints - old_endpoints)
    if len(gains_h) != len(loses_h):
        raise AssertionError("closure rewiring must conserve hydrogen demand")
    nodes = {first_h + index: "H" for index in range(len(gains_h))}
    left = [(offset + source, first_h + index, "-")
            for index, source in enumerate(loses_h)]
    right = [(offset + target, first_h + index, "-")
             for index, target in enumerate(gains_h)]
    return nodes, left, right


def _transition_rule(name: str, attacker: int, old: int, new: int,
                     guard: int | None) -> mf.Rule:
    attacker_offset, target_offset, probe_offset = 0, 40, 80
    attacker_graph = site_heavy(attacker, attacker_offset)
    target_left = site_heavy(old, target_offset)
    target_right = site_heavy(new, target_offset)
    first_connector = _connector(
        attacker_offset + P, target_offset + P, 120
    )
    shared = _merge_overlapping(attacker_graph, target_left, first_connector)
    shared_right = _merge_overlapping(attacker_graph, target_right, first_connector)
    hydrogen_nodes, hydrogen_left, hydrogen_right = _transition_hydrogens(
        old, new, target_offset, 130
    )
    left = make_graph(
        {**node_labels(shared), **hydrogen_nodes},
        (*edge_rows(shared), *hydrogen_left),
    )
    right = make_graph(
        {**node_labels(shared_right), **hydrogen_nodes},
        (*edge_rows(shared_right), *hydrogen_right),
    )
    no_edges = ()
    if guard is not None:
        probe = site_heavy(None, probe_offset)
        second_connector = _connector(
            target_offset + P, probe_offset + P, 140
        )
        left = _merge_overlapping(left, probe, second_connector)
        right = _merge_overlapping(right, probe, second_connector)
        edge = STATE_EDGES[guard]
        no_edges = ((probe_offset + edge[0], probe_offset + edge[1]),)
    return make_rule(name, left, right, no_edges)


def _swap_rule(name: str, before: tuple[int, int], after: tuple[int, int]):
    offsets = (0, 40)
    left_sites = [site_heavy(state, offset) for state, offset in zip(before, offsets)]
    right_sites = [site_heavy(state, offset) for state, offset in zip(after, offsets)]
    connector = _connector(P, 40 + P, 120)
    left = _merge_overlapping(*left_sites, connector)
    right = _merge_overlapping(*right_sites, connector)
    hydrogen_nodes = {}
    hydrogen_left, hydrogen_right = [], []
    next_h = 130
    for old, new, offset in zip(before, after, offsets):
        nodes, old_edges, new_edges = _transition_hydrogens(old, new, offset, next_h)
        hydrogen_nodes.update(nodes)
        hydrogen_left.extend(old_edges)
        hydrogen_right.extend(new_edges)
        next_h += len(nodes)
    return make_rule(
        name,
        make_graph({**node_labels(left), **hydrogen_nodes}, (*edge_rows(left), *hydrogen_left)),
        make_graph({**node_labels(right), **hydrogen_nodes}, (*edge_rows(right), *hydrogen_right)),
    )


def _growth_precursor(state: int, offset=100) -> mf.Graph:
    """Reactive capped repeat used explicitly on the left of an insertion."""
    site = site_heavy(state, offset)
    incoming_o, sulfur, outgoing_o = offset + 10, offset + 11, offset + 12
    first_h, second_h = offset + 13, offset + 14
    return make_graph(
        {
            **node_labels(site),
            incoming_o: "O", sulfur: "S", outgoing_o: "O",
            first_h: "H", second_h: "H",
        },
        (
            *edge_rows(site),
            (incoming_o, offset + P, "-"),
            (offset + P, sulfur, "-"),
            (sulfur, outgoing_o, "-"),
            (incoming_o, first_h, "-"),
            (outgoing_o, second_h, "-"),
        ),
    )


def _growth_rule(name: str, state: int) -> mf.Rule:
    """Net condensation inserting a state-copy repeat into one backbone link."""
    anchor_offset, next_offset = 0, 40
    connector_s, connector_o = 80, 81
    precursor_offset = 100
    ring_context = _merge_overlapping(
        site_heavy(state, anchor_offset),
        site_heavy(None, next_offset),
        _connector(anchor_offset + P, next_offset + P, connector_s),
    )
    precursor = _growth_precursor(state, precursor_offset)
    left = union(ring_context, precursor)

    incoming_o = precursor_offset + 10
    precursor_p = precursor_offset + P
    outgoing_o = precursor_offset + 12
    transferred_h = precursor_offset + 14
    deleted = {
        (min(connector_o, next_offset + P), max(connector_o, next_offset + P), "-"),
        (min(incoming_o, precursor_p), max(incoming_o, precursor_p), "-"),
        (min(outgoing_o, transferred_h), max(outgoing_o, transferred_h), "-"),
    }
    right_edges = [edge for edge in edge_rows(left) if edge not in deleted]
    right_edges.extend(
        (
            (connector_o, precursor_p, "-"),
            (outgoing_o, next_offset + P, "-"),
            (incoming_o, transferred_h, "-"),
        )
    )
    return make_rule(name, left, make_graph(node_labels(left), right_edges))


@lru_cache(None)
def rules() -> dict[str, mf.Rule]:
    result = {
        name: _transition_rule(name, *RULE_TRANSITIONS[name])
        for name in CORE_REACTIONS
    }
    result.update(
        {name: _swap_rule(name, *SWAP_TRANSITIONS[name]) for name in SWAP_REACTIONS}
    )
    result.update(
        {name: _growth_rule(name, GROW_STATES[name]) for name in GROW_REACTIONS}
    )
    return result


def _hydrogen_multiplicity(old: int, new: int) -> int:
    heavy = site_heavy(old)
    sums = bond_order_sums(heavy)
    old_endpoints = set(STATE_EDGES[old][:2])
    new_endpoints = set(STATE_EDGES[new][:2])
    sources = new_endpoints - old_endpoints
    result = 1
    for node in sources:
        result *= VALENCE[node_labels(heavy)[node]] - sums[node]
    return result


@lru_cache(None)
def rate_factors() -> dict[str, Fraction]:
    result = {
        name: Fraction(1, _hydrogen_multiplicity(old, new))
        for name, (_attacker, old, new, _guard) in RULE_TRANSITIONS.items()
    }
    for name, (before, after) in SWAP_TRANSITIONS.items():
        multiplicity = 1
        for old, new in zip(before, after):
            multiplicity *= _hydrogen_multiplicity(old, new)
        result[name] = Fraction(1, multiplicity)
    result.update({name: Fraction(1) for name in GROW_REACTIONS})
    return result


def atom_balance(rule: mf.Rule) -> bool:
    return formula(rule_graph(rule, "left")) == formula(rule_graph(rule, "right"))


def local_bond_order_balance(rule: mf.Rule) -> bool:
    return (node_labels(rule_graph(rule, "left")) == node_labels(rule_graph(rule, "right"))
            and bond_order_sums(rule_graph(rule, "left")) == bond_order_sums(rule_graph(rule, "right")))


def transition_events(states: tuple[int, ...], reaction: str):
    """Exact normalized word events for one reaction channel on a macrocycle."""
    if reaction in GROW_STATES:
        copied = GROW_STATES[reaction]
        states = tuple(states)
        for index, state in enumerate(states):
            if state == copied:
                yield states[:index + 1] + (copied,) + states[index + 1:]
        return
    if reaction in SWAP_TRANSITIONS:
        before, after = SWAP_TRANSITIONS[reaction]
        states = tuple(states)
        for index in range(len(states)):
            if (states[index], states[(index + 1) % len(states)]) != before:
                continue
            target = list(states)
            target[index], target[(index + 1) % len(states)] = after
            yield tuple(target)
        return
    attacker, old, new, guard = RULE_TRANSITIONS[reaction]
    states = tuple(states)
    for index in range(len(states)):
        if states[index] != attacker or states[(index + 1) % len(states)] != old:
            continue
        if guard is not None and states[(index + 2) % len(states)] == guard:
            continue
        target = list(states)
        target[(index + 1) % len(states)] = new
        yield tuple(target)


def word_count(states: tuple[int, ...], word: tuple[int, ...]) -> int:
    states, word = tuple(states), tuple(word)
    return sum(
        all(states[(start + offset) % len(states)] == value
            for offset, value in enumerate(word))
        for start in range(len(states))
    )


def one_site_rhs(state, r: float, g: float):
    """R1 fractions; symmetric swaps and state-copy growth cancel identically."""
    u0, u1, u2 = state
    j01 = r * u0 * u1 + g * u0 * u1 * (1 - u1)
    j12 = r * u1 * u2 + g * u1 * u2 * (1 - u2)
    j20 = r * u2 * u0 + g * u2 * u0 * (1 - u0)
    return j20 - j01, j01 - j12, j12 - j20


def one_site_eigenvalues(r: float, g: float):
    """Physical eigenpair of R1 at u0=u1=u2=1/3."""
    real = -g / 6
    imaginary = r / math.sqrt(3) + g * math.sqrt(3) / 6
    return complex(real, imaginary), complex(real, -imaginary)
