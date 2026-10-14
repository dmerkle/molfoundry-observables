"""Partition enabled source matches by their complete observable jumps."""
from __future__ import annotations
from mfobs.core.graph import edge_rows, node_labels
from mfobs.core.rules import added_edges, preserves_vertices, removed_edges, rule_graph
from collections import defaultdict
from dataclasses import dataclass
import sympy as sp
from mfobs.core.conditions import UnsupportedCondition
from mfobs.models.chlorination.counting import Z, RING, ZERO, ONE, _symbolic_matches, coordinate_boolpoly, certify_locality


@dataclass(frozen=True)
class Event:
    identifier: str
    activity: sp.Expr
    delta: tuple[int, ...]
    certificate: str


def _match_jumps(rule, coordinates):
    """Read the source bond exchange and transport every coordinate per match."""
    def edges(changed):
        return {(frozenset((u, v)), label) for u, v, label in changed}
    if not preserves_vertices(rule):
        raise UnsupportedCondition("Expected a vertex-preserving single-site bond exchange")
    nodes = node_labels(rule_graph(rule, "context"))
    if sorted(nodes.values()) != ["C", "Cl", "H", "H"]:
        raise UnsupportedCondition("Expected a vertex-preserving single-site bond exchange")
    carbon = next(n for n, label in nodes.items() if label == "C")
    chlorine = next(n for n, label in nodes.items() if label == "Cl")
    hydrogens = [n for n, label in nodes.items() if label == "H"]
    pairs = [(h, other) for h in hydrogens for other in hydrogens if h != other
             if edges(removed_edges(rule)) == {(frozenset((carbon, h)), "-"),
                                              (frozenset((other, chlorine)), "-")}]
    if (not pairs or edge_rows(rule_graph(rule, "context")) or edges(added_edges(rule)) != {
            (frozenset((carbon, chlorine)), "-"), (frozenset(hydrogens), "-")}):
        raise UnsupportedCondition("Unsupported source bond exchange")
    for mapping, enabled in _symbolic_matches(rule, 1):
        site = mapping[carbon] - 1
        if site not in range(4):
            raise UnsupportedCondition("Reaction center is outside the neutral ring")
        differences = []
        for before in coordinates.values():
            after = RING.polynomial({powers: c for powers, c in before.terms if not powers[site]})
            differences.append(RING.subtract(after, before))
        yield enabled, differences


def _values(difference):
    # Safe integer bounds from polynomial coefficients, not state evaluation.
    lower = sum(c for _, c in difference.terms if c < 0)
    upper = sum(c for _, c in difference.terms if c > 0)
    if lower.denominator != 1 or upper.denominator != 1:
        raise UnsupportedCondition("Raw count jumps must have integral coefficients")
    return range(int(lower), int(upper) + 1)


def _selector(difference, value):
    if value not in _values(difference):
        return ZERO
    indicator = ONE
    for other in _values(difference):
        if other != value:
            indicator = RING.scale(RING.multiply(indicator, RING.subtract(difference, RING.constant(other))),
                                   sp.Rational(1, value-other))
    return indicator


def _derive_events(rule, reducer, names):
    """Partition actual enabled matches by their entire discovered jump vector."""
    coordinates = {n: reducer.coordinates[n] for n in names}
    counts = defaultdict(lambda: ZERO)
    for enabled, differences in _match_jumps(rule, coordinates):
        cells = {(): enabled}
        for difference in differences:
            cells = {prefix + (value,): selected
                     for prefix, indicator in cells.items()
                     for value in _values(difference)
                     if not (selected := RING.multiply(indicator, _selector(difference, value))).is_zero}
        for delta, indicator in cells.items():
            counts[delta] = RING.add(counts[delta], indicator)
    events = []
    for index, (delta, count) in enumerate(sorted(counts.items()), 1):
        # Support is already fixed: fail closed if joint rates need more counts.
        local = reducer.express(count, discover=False)
        if not local.free_symbols <= {sp.Symbol(n) for n in names}:
            raise UnsupportedCondition("Joint activities require coordinates outside the discovered support")
        events.append(Event(f"event-{index}", sp.expand((reducer.chlorine_total - Z) * local),
                            delta, "derived from symbolic source matches and complete jumps"))
    return tuple(events)


def _certify_events(rule, events, reducer, names):
    """Check disjointness, exhaustiveness and joint rates, not just moments."""
    coordinates = {n: reducer.coordinates[n] for n in names}
    if len({event.delta for event in events}) != len(events):
        raise AssertionError("Event partition has duplicate jump classes")
    counts = [ZERO for _ in events]
    match_count = 0
    for enabled, differences in _match_jumps(rule, coordinates):
        indicators = []
        for event in events:
            if len(event.delta) != len(differences):
                raise AssertionError("Wrong jump dimension")
            indicator = enabled
            for difference, value in zip(differences, event.delta):
                indicator = RING.multiply(indicator, _selector(difference, value))
            if RING.multiply(indicator, indicator) != indicator:
                raise AssertionError("Event selector is not Boolean")
            if any(not RING.multiply(indicator, previous).is_zero for previous in indicators):
                raise AssertionError("Event selectors overlap")
            indicators.append(indicator)
        if RING.add(*indicators) != enabled:
            raise AssertionError("Event partition misses an enabled match")
        counts = [RING.add(count, indicator) for count, indicator in zip(counts, indicators)]
        match_count += 1
    reagent = coordinate_boolpoly(reducer.chlorine_total - Z, coordinates)
    for event, count in zip(events, counts):
        if coordinate_boolpoly(event.activity, coordinates) != RING.multiply(reagent, count):
            raise AssertionError(f"Wrong joint activity for {event.identifier}")
    return {"status": "VERIFIED_SYMBOLIC_JOINT_PARTITION",
            "construction": "symbolic source matches partitioned by complete jump vectors",
            "verification": "literal rule matches, condition AST and all coordinate jumps",
            "disjoint": True, "exhaustive": True,
            "symbolic_matches": match_count, "state_space_enumeration": False}


def _certify_growth(rule, patterns):
    """Certify zero observable jumps by disjointness of changed atom sites."""
    if not preserves_vertices(rule):
        raise UnsupportedCondition("Growth locality requires atom preservation")
    labels = node_labels(rule_graph(rule, "left"))
    left = set(edge_rows(rule_graph(rule, "left")))
    changed = removed_edges(rule) ^ added_edges(rule)
    outside = {u for u, label in labels.items() if label == "C" and
               sum(u in (a, b) and bond == "-" and labels[b if a == u else a] == "H"
                   for a, b, bond in left) >= 2}
    outside |= {u for u, label in labels.items() if label == "H" and
                any(u in (a, b) and bond == "-" and (b if a == u else a) in outside
                    for a, b, bond in left)}
    if any(bond != "-" or u not in outside or v not in outside for u, v, bond in changed):
        raise UnsupportedCondition("Growth may change the ring/HCl core")
    for pattern in patterns.values():
        certify_locality(pattern.rule)
    return {"status": "STRUCTURAL_ZERO_JUMP",
            "changed_vertices_outside_ring_and_HCl": sorted({u for a, b, _ in changed for u in (a, b)}),
            "argument": "Every changed carbon has at least two positive H neighbours, so cannot be "
                        "a ring carbon; changed hydrogens belong to those carbons. Ring bonds, "
                        "substituents and HCl are preserved. All discovered positive patterns are core-local.",
            "all_tail_lengths": True, "state_space_enumeration": False}
