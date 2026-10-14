"""Bond-only chemical graphs for a Brusselator coupled to unbounded ROP.

This is a deliberately hypothetical net-reaction benchmark.  The Brusselator
species are represented by explicit-H constitutional isomers, so their identity
is carried entirely by ordinary covalent connectivity.  No ``X``/``Y`` atom
labels, state edges, isotopes, or charges are used.

The four core channels preserve every atom.  Chemostatting A and B and removing
D and E is an open-reactor boundary assumption made only when projecting the
graph rules to the deterministic Brusselator ODE.  The separate monomer inlet is
an explicit atom-creating boundary operation, as in the existing ROP example.
These rules are not proposed as elementary hydrocarbon mechanisms.
"""
from __future__ import annotations
from mfobs.core.graph import make_graph, node_labels
from mfobs.core.rules import make_rule, rule_graph
import molfoundry as mf

from fractions import Fraction
from functools import lru_cache


from mfobs.core.molecules import automorphism_count, formula
from mfobs.models.molecular_rop import inputs as rop
from mfobs.core.gml import bond_order_sums, disjoint_union, shift, union
from mfobs.core.molecules import complete_with_h
from mfobs.state_space import StateSpace, covalent_theory


def state_family():
    """Declare Ω for this model, independently of its selected queries."""
    return StateSpace("brusselator-growth", {
        "molecules": ["A", "B", "D", "E", "X", "Y"],
        "dynamic_coordinates": {"X": "X", "Y": "Y"},
        "persistent_components": ["A", "B"], "discarded_components": ["D", "E"],
        "polymer_alphabet": ["B"],
        "path": "one finite caprolactone chain and arbitrary finite monomer pool",
        "boundary": "A/B chemostats and removal of D/E in the open-reactor projection",
    }, covalent_theory())


SOURCE = "source"
AUTOCATALYSIS = "autocatalysis"
CONVERSION = "conversion"
SINK = "sink"
MONOMER_INLET = "monomer_inlet"
GROWTH = "growth"
REACTION_ORDER = (
    SOURCE,
    AUTOCATALYSIS,
    CONVERSION,
    SINK,
    MONOMER_INLET,
    GROWTH,
)
RATE_NAMES = {
    SOURCE: "k1",
    AUTOCATALYSIS: "k2",
    CONVERSION: "k3",
    SINK: "k4",
    MONOMER_INLET: "beta",
    GROWTH: "k_g",
}
RATE_CHANNELS = {
    SOURCE: "k1*a*V",
    AUTOCATALYSIS: "k2/V^2",
    CONVERSION: "k3*b",
    SINK: "k4",
    MONOMER_INLET: "beta",
    GROWTH: "k_g",
}


def _hydrocarbon(size: int, edges) -> mf.Graph:
    graph = complete_with_h(make_graph({u: "C" for u in range(size)}, edges))
    state_family().validate_graph(graph)
    return graph


@lru_cache(None)
def species(name: str) -> mf.Graph:
    """Return one explicit-H molecular graph for a named formal species.

    A, X, Y, and E are four C4H8 isomers.  B and D are the two graph-distinct
    C3H6 isomers.  The assignments are encoding choices, not a mechanistic
    statement about real isomerization kinetics.
    """
    builders = {
        # cyclobutane
        "A": lambda: _hydrocarbon(
            4, ((0, 1, "-"), (1, 2, "-"), (2, 3, "-"), (3, 0, "-"))
        ),
        # 1-butene
        "X": lambda: _hydrocarbon(
            4, ((0, 1, "="), (1, 2, "-"), (2, 3, "-"))
        ),
        # methylcyclopropane
        "Y": lambda: _hydrocarbon(
            4, ((0, 1, "-"), (1, 2, "-"), (2, 0, "-"), (2, 3, "-"))
        ),
        # isobutene
        "E": lambda: _hydrocarbon(
            4, ((0, 1, "="), (0, 2, "-"), (0, 3, "-"))
        ),
        # cyclopropane
        "B": lambda: _hydrocarbon(
            3, ((0, 1, "-"), (1, 2, "-"), (2, 0, "-"))
        ),
        # propene
        "D": lambda: _hydrocarbon(
            3, ((0, 1, "="), (1, 2, "-"))
        ),
    }
    try:
        return builders[name]()
    except KeyError as exc:
        raise ValueError(f"Unknown Brusselator species {name!r}") from exc


@lru_cache(None)
def boundary_anchor() -> mf.Graph:
    """One unchanged water molecule localizes the monomer boundary rule."""
    return make_graph(
        {0: "O", 1: "H", 2: "H"}, ((0, 1, "-"), (0, 2, "-"))
    )


def _mapped_union(*names: str) -> mf.Graph:
    return disjoint_union([species(name) for name in names])


@lru_cache(None)
def core_rules() -> dict[str, mf.Rule]:
    """Four atom-preserving bond-rewrite channels of the Brusselator."""
    return {
        SOURCE: make_rule(
            "brusselator-A-to-X", species("A"), species("X")
        ),
        AUTOCATALYSIS: make_rule(
            "brusselator-2X-plus-Y-to-3X",
            _mapped_union("X", "X", "Y"),
            _mapped_union("X", "X", "X"),
        ),
        CONVERSION: make_rule(
            "brusselator-B-plus-X-to-Y-plus-D",
            _mapped_union("B", "X"),
            _mapped_union("D", "Y"),
        ),
        SINK: make_rule(
            "brusselator-X-to-E", species("X"), species("E")
        ),
    }


@lru_cache(None)
def monomer_inlet_rule() -> mf.Rule:
    anchor = boundary_anchor()
    return make_rule(
        "caprolactone-inlet-at-water-boundary",
        anchor,
        union(anchor, shift(rop.monomer("B"), 3)),
    )


@lru_cache(None)
def growth_rule() -> mf.Rule:
    """X-catalyzed caprolactone addition; X is unchanged graph context."""
    left, right = rop.propagation_sides("B")
    catalyst = shift(species("X"), 10)
    return make_rule(
        "X-catalyzed-caprolactone-net-propagation",
        union(left, catalyst),
        union(right, catalyst),
    )


@lru_cache(None)
def rules() -> dict[str, mf.Rule]:
    result = {**core_rules(), MONOMER_INLET: monomer_inlet_rule(), GROWTH: growth_rule()}
    omega = state_family()
    for rule in result.values():
        omega.validate_rule(rule)
    return result


@lru_cache(None)
def rate_factors() -> dict[str, Fraction]:
    """Normalize raw injective maps to the displayed mass-action channels.

    For 2X+Y the factor intentionally retains the ordered falling factorial
    (N_X)_2 N_Y, whose density limit is x^2 y.  Dividing by the whole left-side
    automorphism group would instead introduce an unwanted factor 1/2.
    """
    aut = {name: automorphism_count(species(name)) for name in "AXYEBD"}
    return {
        SOURCE: Fraction(1, aut["A"]),
        AUTOCATALYSIS: Fraction(1, aut["X"] ** 2 * aut["Y"]),
        CONVERSION: Fraction(1, aut["B"] * aut["X"]),
        SINK: Fraction(1, aut["X"]),
        MONOMER_INLET: Fraction(1, automorphism_count(boundary_anchor())),
        GROWTH: Fraction(1, aut["X"]),
    }


def atom_balance(rule: mf.Rule) -> bool:
    return formula(rule_graph(rule, "left")) == formula(rule_graph(rule, "right"))


def local_bond_order_balance(rule: mf.Rule) -> bool:
    """Every preserved atom has the same ordinary-covalent valence on both sides."""
    if node_labels(rule_graph(rule, "left")) != node_labels(rule_graph(rule, "right")):
        return False
    return bond_order_sums(rule_graph(rule, "left")) == bond_order_sums(rule_graph(rule, "right"))


def event_rates(parameters, volume: float) -> dict[str, float]:
    """Translate density parameters to the six exported count-rate channels."""
    if volume <= 0:
        raise ValueError("volume must be positive")
    return {
        RATE_CHANNELS[SOURCE]: parameters.get("k1", 1.0) * parameters["a"] * volume,
        RATE_CHANNELS[AUTOCATALYSIS]: parameters.get("k2", 1.0) / volume**2,
        RATE_CHANNELS[CONVERSION]: parameters.get("k3", 1.0) * parameters["b"],
        RATE_CHANNELS[SINK]: parameters.get("k4", 1.0),
        RATE_CHANNELS[MONOMER_INLET]: parameters["beta"],
        RATE_CHANNELS[GROWTH]: parameters["k_g"] / volume,
    }


def deterministic_rhs(_time: float, state, parameters) -> tuple[float, ...]:
    """Held-out macroscopic reference in coordinates (x, y, M, n)."""
    x, y, monomer, length = state
    del length
    a = parameters["a"]
    b = parameters["b"]
    k1 = parameters.get("k1", 1.0)
    k2 = parameters.get("k2", 1.0)
    k3 = parameters.get("k3", 1.0)
    k4 = parameters.get("k4", 1.0)
    beta = parameters["beta"]
    growth = parameters["k_g"] * x * monomer
    return (
        k1 * a + k2 * x * x * y - k3 * b * x - k4 * x,
        k3 * b * x - k2 * x * x * y,
        beta - growth,
        growth,
    )


def core_steady_state(parameters) -> tuple[float, float]:
    """Steady state for the unit-rate Brusselator used by the walkthrough."""
    if any(parameters.get(key, 1.0) != 1.0 for key in ("k1", "k2", "k3", "k4")):
        raise ValueError("Closed-form steady state is documented only for unit core rates")
    return float(parameters["a"]), float(parameters["b"] / parameters["a"])


def hopf_threshold(a: float) -> float:
    return 1.0 + a * a


def core_jacobian_at_steady(parameters) -> tuple[tuple[float, float], ...]:
    if any(parameters.get(key, 1.0) != 1.0 for key in ("k1", "k2", "k3", "k4")):
        raise ValueError("Closed-form Jacobian is documented only for unit core rates")
    a, b = float(parameters["a"]), float(parameters["b"])
    return ((b - 1.0, a * a), (-b, -a * a))
