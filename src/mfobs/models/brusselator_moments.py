"""Exact factorial-moment discovery for the bond-isomer Brusselator core.

The chemical input remains a collection of ordinary element/bond graph rules.
This adapter first reads their X/Y population changes from the literal GML and
then applies the exact reaction-channel generator to factorial count
observables ``(N_X)_p (N_Y)_q``.  Each newly required coordinate is exported as
an identity rule on ``p`` disconnected X molecules and ``q`` disconnected Y
molecules.

This is deliberately a population-language reducer, not a MolFoundry overlap
enumeration of the very large explicit-H patterns.  No mean-field
factorisation is used.  Since the autocatalytic channel raises moment degree,
the shared worklist is normally stopped by an exploration budget and reports
UNKNOWN while retaining every observable found before the cutoff.
"""
from __future__ import annotations
from mfobs.core.graph import make_graph
from mfobs.core.rules import make_rule
import molfoundry as mf

from fractions import Fraction
from math import prod

import sympy as sp


from mfobs.core.gml import rule_side
from mfobs.core.molecules import automorphism_count
from mfobs.models.population_projection import molecular_component_inventory

from mfobs.core.observables import parse_observable
from mfobs.core.semantic_contract import SemanticContract, theorem_preservation
from mfobs.inputs import InputError


NX, NY = sp.symbols("N_X N_Y", integer=True, nonnegative=True)
VOLUME, A_FEED, B_FEED = sp.symbols("V a b", positive=True)

CHANNELS = {
    "source": {"delta": (1, 0), "propensity": A_FEED * VOLUME},
    "autocatalysis": {
        "delta": (1, -1),
        "propensity": sp.ff(NX, 2) * NY / VOLUME**2,
    },
    "conversion": {"delta": (-1, 1), "propensity": B_FEED * NX},
    "sink": {"delta": (-1, 0), "propensity": NX},
}

EXPECTED_STOICHIOMETRY = {
    "source": ({"A": 1}, {"X": 1}),
    "autocatalysis": ({"X": 2, "Y": 1}, {"X": 3}),
    "conversion": ({"B": 1, "X": 1}, {"D": 1, "Y": 1}),
    "sink": ({"X": 1}, {"E": 1}),
}


def _clean(form):
    return {
        key: sp.factor(value)
        for key, value in sorted(form.items())
        if sp.cancel(value) != 0
    }


def factorial_expression(moment: tuple[int, int]):
    """Return the exact falling-factorial polynomial for one coordinate."""
    p, q = moment
    return sp.expand(sp.ff(NX, p) * sp.ff(NY, q))


def ordinary_to_factorial(expression):
    """Expand a polynomial in N_X,N_Y into the unique factorial basis."""
    polynomial = sp.Poly(sp.expand(expression), NX, NY)
    result = {}
    for (x_degree, y_degree), coefficient in polynomial.terms():
        for p in range(x_degree + 1):
            sx = sp.functions.combinatorial.numbers.stirling(
                x_degree, p, kind=2
            )
            if not sx:
                continue
            for q in range(y_degree + 1):
                sy = sp.functions.combinatorial.numbers.stirling(
                    y_degree, q, kind=2
                )
                value = coefficient * sx * sy
                if value:
                    result[(p, q)] = result.get((p, q), sp.S.Zero) + value
    return _clean(result)


def channel_derivative(moment: tuple[int, int], channel: str):
    """Exact generator image of one factorial moment in one rate channel."""
    try:
        record = CHANNELS[channel]
    except KeyError as exc:
        raise InputError(f"Unsupported Brusselator channel {channel!r}") from exc
    delta_x, delta_y = record["delta"]
    observable = factorial_expression(moment)
    jumped = observable.subs(
        {NX: NX + delta_x, NY: NY + delta_y}, simultaneous=True
    )
    return ordinary_to_factorial(record["propensity"] * (jumped - observable))


def moment_name(moment: tuple[int, int]) -> str:
    p, q = moment
    if (p, q) == (0, 0):
        return "ONE"
    x = "" if p == 0 else ("X" if p == 1 else f"X{p}")
    y = "" if q == 0 else ("Y" if q == 1 else f"Y{q}")
    return x + y


def moment_meaning(moment: tuple[int, int]) -> str:
    p, q = moment
    factors = []
    if p:
        factors.append("N_X" if p == 1 else f"(N_X)_{p}")
    if q:
        factors.append("N_Y" if q == 1 else f"(N_Y)_{q}")
    return "1" if not factors else " ".join(factors)


def _graph(record) -> mf.Graph:
    return make_graph(record["nodes"], record["edges"])


class Adapter:
    """Worklist adapter for exact, dynamically discovered factorial moments."""

    semantic_contract_required = True

    def __init__(self, inputs):
        inputs.validate_overlap_policy()
        configuration = inputs.read("reduction")
        if set(configuration) != {"method", "coordinate_family"}:
            raise InputError("Unexpected Brusselator moment reduction fields")
        if configuration["method"] != "exact-factorial-population-generator-v1":
            raise InputError("Wrong Brusselator moment reduction method")
        if configuration["coordinate_family"] != "all-(N_X)_p-(N_Y)_q":
            raise InputError("Wrong Brusselator factorial coordinate family")

        specification = inputs.read("state_language")
        if specification.get("kind") != "open-brusselator-X-Y-population-language":
            raise InputError("Wrong Brusselator population state language")
        finite = specification.get("finite_components", {})
        if not {"A", "B", "D", "E", "X", "Y"} <= set(finite):
            raise InputError("State language is missing a Brusselator component")

        self.inputs = inputs
        self.specification = specification
        self.species = {name: _graph(record) for name, record in finite.items()}
        self.automorphisms = {
            name: automorphism_count(graph) for name, graph in self.species.items()
        }
        self.moments = {}
        self.names = {}
        self.origins = {}
        self.discoveries = []
        self.generator_rows = {}

        actual_stoichiometry = self._read_stoichiometry()
        if actual_stoichiometry != EXPECTED_STOICHIOMETRY:
            raise InputError(
                "Literal GML rule stoichiometry does not match the "
                f"Brusselator channel contract: {actual_stoichiometry}"
            )
        self.stoichiometry = actual_stoichiometry
        expected_per_match = {
            name: Fraction(1, prod(self.automorphisms[species] ** count
                                  for species, count in left.items()))
            for name, (left, _right) in actual_stoichiometry.items()
        }
        supplied_per_match = {
            reaction["id"]: Fraction(reaction["per_match"])
            for reaction in inputs.data["rules"]
        }
        if supplied_per_match != expected_per_match:
            raise InputError(
                "GML per-match normalisation does not recover the declared "
                "population propensities"
            )
        self.per_match_normalisation = {
            name: str(value) for name, value in expected_per_match.items()
        }
        scope = "all nonnegative integer X/Y populations at fixed V>0"
        preservation = tuple(
            theorem_preservation(
                reaction["id"],
                scope,
                {
                    "literal_gml_stoichiometry": [
                        self.stoichiometry[reaction["id"]][0],
                        self.stoichiometry[reaction["id"]][1],
                    ],
                    "boundary_argument": (
                        "Every negative population change is multiplied by its "
                        "required falling-factorial reactant count, hence its "
                        "propensity vanishes at the boundary."
                    ),
                },
            )
            for reaction in inputs.data["rules"]
        )
        self.semantic_contract = SemanticContract(
            domain=specification["kind"],
            reducer="exact symbolic finite differences in the factorial basis",
            exact_reduction_scope=(
                "Identity patterns made of complete, disconnected X/Y molecular "
                "components; coefficients lie in Q(a,b,V)."
            ),
            preservation=preservation,
        )

    def _component_inventory(self, graph: mf.Graph):
        try:
            return molecular_component_inventory(graph, self.species)
        except ValueError as exc:
            raise InputError(str(exc)) from exc

    def _read_stoichiometry(self):
        result = {}
        rules = self.inputs.engine_rules()
        if set(rules) != set(CHANNELS):
            raise InputError(f"Expected exactly the core channels {sorted(CHANNELS)}")
        for name, rule in rules.items():
            result[name] = (
                self._component_inventory(rule_side(rule, "left")),
                self._component_inventory(rule_side(rule, "right")),
            )
        return result

    def _moment_from_rule(self, rule) -> tuple[int, int]:
        parsed = parse_observable(rule)
        if not parsed.jump_identity or parsed.no_edges or parsed.generic_conditions:
            raise InputError("Moment seeds must be unconditional identity rules")
        graph = make_graph(parsed.labels, parsed.edges)
        inventory = self._component_inventory(graph)
        if set(inventory) - {"X", "Y"}:
            raise InputError("Moment seeds may contain only complete X/Y molecules")
        return inventory.get("X", 0), inventory.get("Y", 0)

    def _identity_rule(self, moment: tuple[int, int]):
        from mfobs.core.gml import disjoint_union

        p, q = moment
        graph = disjoint_union(
            [self.species["X"] for _ in range(p)]
            + [self.species["Y"] for _ in range(q)]
        )
        identity = make_rule(moment_name(moment), graph, graph)
        return identity

    def _normalisation(self, moment: tuple[int, int]) -> Fraction:
        p, q = moment
        return Fraction(
            self.automorphisms["X"] ** p * self.automorphisms["Y"] ** q,
            1,
        )

    def _add_moment(self, moment, state, *, origin, supplied_rule=None):
        if moment in self.names:
            return self.names[moment]
        name = moment_name(moment)
        if name in state.patterns:
            raise InputError(f"Moment name collision for {moment}: {name}")
        rule = supplied_rule if supplied_rule is not None else self._identity_rule(moment)
        state.add(name, rule, moment_meaning(moment))
        self.moments[name] = moment
        self.names[moment] = name
        self.origins[name] = origin
        return name

    def initialize(self, state, composer):
        observables = self.inputs.engine_observables()
        for name in self.inputs.data["initial_observables"]:
            rule = observables[name]
            moment = self._moment_from_rule(rule)
            if moment == (0, 0):
                raise InputError("The initial moment must be nonconstant")
            expected = moment_name(moment)
            if name != expected:
                raise InputError(f"Seed {name!r} must use semantic name {expected!r}")
            self._add_moment(
                moment,
                state,
                supplied_rule=rule,
                origin={
                    "kind": "input-seed",
                    "file": next(
                        item["file"]
                        for item in self.inputs.data["observables"]
                        if item["id"] == name
                    ),
                },
            )

    def expand(self, name, state, composer):
        moment = self.moments[name]
        state.equations[name] = {}
        self.generator_rows[name] = {}
        for reaction in self.inputs.data["rules"]:
            channel = reaction["id"]
            derivative = channel_derivative(moment, channel)
            row = {}
            source_scale = self._normalisation(moment)
            for target_moment, coefficient in sorted(
                derivative.items(), key=lambda item: (sum(item[0]), item[0])
            ):
                is_new = target_moment not in self.names
                target = self._add_moment(
                    target_moment,
                    state,
                    origin={
                        "kind": "exact-generator-residual",
                        "source": name,
                        "reaction": channel,
                        "factorial_moment": list(target_moment),
                    },
                )
                target_scale = self._normalisation(target_moment)
                raw_coefficient = coefficient * sp.Rational(
                    source_scale.numerator * target_scale.denominator,
                    source_scale.denominator * target_scale.numerator,
                )
                row[target] = str(sp.factor(raw_coefficient))
                if is_new:
                    self.discoveries.append(
                        {
                            "source": name,
                            "reaction": channel,
                            "added": target,
                            "moment": list(target_moment),
                            "factorial_coefficient": str(sp.factor(coefficient)),
                            "raw_graph_coefficient": row[target],
                        }
                    )
            state.equations[name][channel] = row
            self.generator_rows[name][channel] = {
                "source_moment": list(moment),
                "terms": {
                    moment_name(target): str(value)
                    for target, value in derivative.items()
                },
                "method": "propensity * (f(N+delta)-f(N)), exact factorial conversion",
            }

    def finalize(self, state, composer, complete):
        initial = self.inputs.data["initial_state"]
        if set(initial) != {"X", "Y", "V"}:
            raise InputError("initial_state must contain exactly X, Y and V")
        if any(type(initial[key]) is not int or initial[key] < 0 for key in ("X", "Y")):
            raise InputError("Initial X/Y populations must be nonnegative integers")
        if type(initial["V"]) is not int or initial["V"] <= 0:
            raise InputError("Initial V must be a positive integer")

        initial_values = {}
        normalisation = {}
        coordinate_records = {}
        for name, (p, q) in self.moments.items():
            scale = self._normalisation((p, q))
            initial_values[name] = str(
                scale * sp.ff(initial["X"], p) * sp.ff(initial["Y"], q)
            )
            normalisation[name] = {
                "scale": str(scale),
                "meaning": moment_meaning((p, q)),
            }
            coordinate_records[name] = {
                "p": p,
                "q": q,
                "degree": p + q,
                "factorial_moment": moment_meaning((p, q)),
                "raw_graph_count_scale": str(scale),
            }

        return {
            "status": (
                "CLOSED_EXACT_FACTORIAL_MOMENT_SYSTEM" if complete else "UNKNOWN"
            ),
            "proof_status": (
                "Exact reaction-channel finite differences and exact Stirling "
                "conversion to graph-represented factorial moments; no moment "
                "factorisation. Literal GML fixes channel stoichiometry."
            ),
            "coefficient_parameters": ["V", "a", "b"],
            "parameter_domain": {
                "V": "positive system-size parameter",
                "a": "nonnegative A-feed concentration",
                "b": "nonnegative B-feed concentration",
            },
            "population_stoichiometry_from_gml": self.stoichiometry,
            "per_match_normalisation_verified": self.per_match_normalisation,
            "moment_coordinates": coordinate_records,
            "basis_origins": self.origins,
            "discovery_decisions": self.discoveries,
            "generator_derivations": self.generator_rows,
            "normalisation": normalisation,
            "initial_values": initial_values,
            "synthesis_state_evaluations": 0,
            "synthesis_graph_compositions": 0,
            "composition_note": (
                "Large explicit-H overlap composition is replaced by the exact "
                "population generator after its stoichiometry is read from the "
                "literal chemical GML rules."
            ),
            "basis_kind": (
                "complete-molecule identity patterns representing independent "
                "X/Y falling-factorial moments"
            ),
            "limitations": [
                "An expansion cutoff is UNKNOWN, not a finite exact closure.",
                "The adapter supports this two-population mass-action state language, not arbitrary graph grammars.",
                "Chemostatted A/B and removed D/E are open-reactor assumptions.",
                "The deterministic nonlinear ODE shown in the report is a large-volume limit, not the exported exact expectation system.",
                "The bond-isomer rules are a hypothetical net-reaction benchmark, not an elementary hydrocarbon mechanism.",
            ],
        }
