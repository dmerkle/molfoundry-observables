"""Ring-8 finite-state oracle and shortest-word normal forms (not chemistry)."""
from __future__ import annotations
from mfobs.core.graph import make_graph, observable
from mfobs.core.engine import mf
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
import sympy as sp
from mfobs.core.rules import rule_graph, rule_vertex_ids
from mfobs.core.conditions import identity_rule
from mfobs.core.gml import graph_gml
BitState=tuple[int,...]
def ring_host_gml(state: tuple[int, ...]) -> str:
    return graph_gml(slot_ring_graph(state))


def slot_ring_graph(state: tuple[int, ...]) -> mf.Graph:
    return _slot_graph(state, close_cycle=True)


def _slot_graph(state: tuple[int, ...], *, close_cycle: bool) -> mf.Graph:
    nodes: dict[int, str] = {}
    edges: list[tuple[int, int, str]] = []
    n = len(state)
    for site, occupied in enumerate(state):
        center, port, ligand = 3 * site, 3 * site + 1, 3 * site + 2
        nodes.update({center: "C", port: "P", ligand: "L"})
        edges.extend([(center, port, "cp"), (center, ligand, "cl")])
        if occupied:
            edges.append((port, ligand, "bound"))
    for site in range(n if close_cycle else n - 1):
        edges.append((3 * site, 3 * ((site + 1) % n), "next"))
    return make_graph(nodes, edges)

def dihedral_images(state: BitState) -> set[BitState]:
    n = len(state)
    reflected = state[::-1]
    return {
        *(state[i:] + state[:i] for i in range(n)),
        *(reflected[i:] + reflected[:i] for i in range(n)),
    }

def canonical_orbit(state: BitState) -> BitState:
    return min(dihedral_images(state))

@dataclass(frozen=True)
class RingRates:
    basal_bind: int = 1
    cooperative_bind: int = 1
    unbind: int = 1

def transitions(state: BitState, rates: RingRates = RingRates()):
    """Yield flips and rates induced by add-edge NAC/delete-edge rules.

    An empty slot binds basally and once for each occupied neighbour.  These
    are distinct rule matches, exactly as in a stochastic graph grammar.
    """

    n = len(state)
    for i, occupied in enumerate(state):
        target = list(state)
        target[i] = 1 - occupied
        if occupied:
            rate = rates.unbind
        else:
            occupied_neighbours = state[(i - 1) % n] + state[(i + 1) % n]
            rate = rates.basal_bind + rates.cooperative_bind * occupied_neighbours
        yield tuple(target), rate

def orbit_generator(n: int, rates: RingRates = RingRates()):
    """Construct and independently verify the exact dihedral lumping."""

    states = list(product((0, 1), repeat=n))
    orbits = sorted({canonical_orbit(state) for state in states})
    index = {orbit: i for i, orbit in enumerate(orbits)}
    generator = sp.zeros(len(orbits), len(orbits))

    for row, representative in enumerate(orbits):
        reference: dict[int, int] | None = None
        for state in dihedral_images(representative):
            aggregated: defaultdict[int, int] = defaultdict(int)
            total = 0
            for target, rate in transitions(state, rates):
                aggregated[index[canonical_orbit(target)]] += rate
                total += rate
            if reference is None:
                reference = dict(aggregated)
            elif dict(aggregated) != reference:
                raise AssertionError("the proposed orbit partition is not lumpable")
        assert reference is not None
        for column, rate in reference.items():
            generator[row, column] += rate
        generator[row, row] -= sum(rate for _, rate in transitions(representative, rates))
    return orbits, generator
def word_observable(word: tuple[int, ...], *, close_cycle: bool = False):
    """Identity observable for an occupied/empty word on a ring path."""
    if not word:
        raise ValueError("a motif word must be non-empty")
    if close_cycle and len(word) < 3:
        raise ValueError("cycle word must contain at least three slots")
    graph = _slot_graph(word, close_cycle=close_cycle)
    absent = tuple((3 * site + 1, 3 * site + 2)
                   for site, occupied in enumerate(word) if not occupied)
    name = ("cycle" if close_cycle else "word") + "-" + "".join(map(str, word))
    return identity_rule(graph, name, absent)


def word_values(
    word: tuple[int, ...], orbits: list[tuple[int, ...]]
) -> tuple[int, ...]:
    """Count path embeddings by cyclic starts and orientations."""

    result = []
    for state in orbits:
        n = len(state)
        count = 0
        for start in range(n):
            for direction in (-1, 1):
                if all(
                    state[(start + direction * offset) % n] == bit
                    for offset, bit in enumerate(word)
                ):
                    count += 1
        # A one-vertex path has no orientation; the graph matcher counts each
        # slot once rather than twice.
        if len(word) == 1:
            count //= 2
        result.append(count)
    return tuple(result)

@dataclass(frozen=True)
class ObservableSpec:
    input_graph: mf.Graph
    no_edges: tuple[tuple[int, int], ...]
    dpo_deleted_nodes: tuple[int, ...]

def observable_spec(rule) -> ObservableSpec:
    """Extract the input pattern and supported conditions from a jump span."""

    if rule.matchCondition is not None:
        raise NotImplementedError("generic nested conditions are not yet evaluated")
    if rule.matchConstraintsAdjacency:
        raise NotImplementedError("adjacency constraints are not yet evaluated")
    if rule.matchConstraintsLabelAny:
        raise NotImplementedError("label-any constraints are not yet evaluated")
    if rule.matchConstraintsShortestPath:
        raise NotImplementedError("shortest-path constraints are not yet evaluated")

    return ObservableSpec(
        rule_graph(rule),
        tuple((c.source, c.target) for c in rule.matchConstraintsNoEdge),
        rule_vertex_ids(rule, "deleted"),
    )

def _evaluate_spec(spec: ObservableSpec, states: tuple[mf.Graph, ...]) -> tuple[int, ...]:
    if spec.dpo_deleted_nodes:
        # The present ring grammar preserves every node, so this branch is not
        # needed.  Silently treating a DPO q:K->I with deleted nodes as an
        # ordinary pattern count would lose its dangling-condition semantics.
        raise NotImplementedError("DPO observables deleting nodes need a gluing check")
    return tuple(
        observable(spec.input_graph, state, spec.no_edges) for state in states
    )

def evaluate_observable(rule, states: list[mf.Graph]) -> tuple[int, ...]:
    return _evaluate_spec(observable_spec(rule), tuple(states))
