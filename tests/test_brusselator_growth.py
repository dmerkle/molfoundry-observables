"""Checks for the bond-isomer Brusselator with unbounded polymer growth."""
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import edge_rows, make_graph, node_labels, node_rows, same_graph
from mfobs.core.rules import rule_graph
from fractions import Fraction

import numpy as np
import pytest

from mfobs.core.conditions import matches
from mfobs.core.graph import apply_match
from mfobs.models import brusselator_growth_chemistry as chemistry
from mfobs.models import brusselator_moments as moments
from mfobs.models.molecular_rop import inputs as rop
from mfobs.core.gml import bond_order_sums, disjoint_union
from mfobs.core.molecules import VALENCE, formula


def test_species_are_explicit_valid_bond_isomers_not_state_labels():
    expected = {
        "A": {"C": 4, "H": 8}, "X": {"C": 4, "H": 8},
        "Y": {"C": 4, "H": 8}, "E": {"C": 4, "H": 8},
        "B": {"C": 3, "H": 6}, "D": {"C": 3, "H": 6},
    }
    for name, expected_formula in expected.items():
        graph = chemistry.species(name)
        assert chemistry.formula(graph) == expected_formula
        assert set(node_labels(graph).values()) == {"C", "H"}
        assert {label for _, _, label in edge_rows(graph)} <= {"-", "="}
        sums = bond_order_sums(graph)
        assert all(sums[u] == VALENCE[label] for u, label in node_rows(graph))


def test_rdkit_accepts_all_full_molecular_components():
    Chem = pytest.importorskip("rdkit.Chem")
    for graph in [chemistry.species(name) for name in "AXYEBD"] + [
        chemistry.boundary_anchor(), rop.monomer("B"), rop.chain(()), rop.chain(("B",))
    ]:
        assert {label for _, _, label in edge_rows(graph)} <= {"-", "="}
        molecule = to_rdkit(graph)
        assert molecule.GetNumAtoms() == graph.numVertices
        assert all(atom.GetNumImplicitHs() == 0 for atom in molecule.GetAtoms())
        # RDKit aromaticizes the benzyl ring after parsing its Kekule -/= input.
        assert all(bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                                         Chem.BondType.AROMATIC)
                   for bond in molecule.GetBonds())


def test_core_rule_sides_preserve_atoms_and_local_valence():
    for rule in chemistry.core_rules().values():
        assert chemistry.atom_balance(rule)
        assert chemistry.local_bond_order_balance(rule)
        assert node_rows(rule_graph(rule, "left")) == node_rows(rule_graph(rule, "right"))
        assert {label for _, _, label in (*edge_rows(rule_graph(rule, "left")), *edge_rows(rule_graph(rule, "right")))} <= {"-", "="}
        identity = dict(node_labels(rule_graph(rule, "left")))
        identity = {node: node for node in identity}
        assert same_graph(apply_match(rule, rule_graph(rule, "left"), identity), rule_graph(rule, "right"))


def test_growth_is_atom_preserving_and_inlet_is_the_only_boundary_creation():
    growth = chemistry.growth_rule()
    assert chemistry.atom_balance(growth)
    assert chemistry.local_bond_order_balance(growth)
    inlet = chemistry.monomer_inlet_rule()
    assert not chemistry.atom_balance(inlet)
    created = set(node_labels(rule_graph(inlet, "right"))) - set(node_labels(rule_graph(inlet, "left")))
    output = rule_graph(inlet, "right")
    monomer = make_graph({u: label for u, label in node_rows(output) if u in created},
                         [(u, v, label) for u, v, label in edge_rows(output) if u in created and v in created])
    assert formula(monomer) == formula(rop.monomer("B"))
    for name, rule in chemistry.rules().items():
        assert chemistry.atom_balance(rule) == (name != chemistry.MONOMER_INLET)


def test_per_match_factors_recover_population_mass_action_counts():
    factors = chemistry.rate_factors()
    source = chemistry.core_rules()[chemistry.SOURCE]
    assert factors[chemistry.SOURCE] * len(matches(rule_graph(source, "left"), rule_graph(source, "left"))) == 1

    conversion = chemistry.core_rules()[chemistry.CONVERSION]
    assert factors[chemistry.CONVERSION] * len(matches(rule_graph(conversion, "left"), rule_graph(conversion, "left"))) == 1

    autocatalysis = chemistry.core_rules()[chemistry.AUTOCATALYSIS]
    # The two indistinguishable X components deliberately give (2)_2 = 2.
    assert factors[chemistry.AUTOCATALYSIS] * len(
        matches(rule_graph(autocatalysis, "left"), rule_graph(autocatalysis, "left"))
    ) == 2

    inlet = chemistry.monomer_inlet_rule()
    assert factors[chemistry.MONOMER_INLET] * len(matches(rule_graph(inlet, "left"), rule_graph(inlet, "left"))) == 1

    growth = chemistry.growth_rule()
    host = disjoint_union([rop.chain(()), rop.monomer("B"), chemistry.species("X")])
    raw_growth = len(matches(growth, host))
    assert factors[chemistry.GROWTH] * raw_growth == 1
    assert factors[chemistry.GROWTH] == Fraction(1, 24)


@pytest.mark.parametrize("b, expected_sign", [(1.8, -1), (2.2, 1)])
def test_hopf_threshold_and_core_linearization(b, expected_sign):
    parameters = {"a": 1.0, "b": b, "beta": 0.3, "k_g": 0.15}
    steady = (*chemistry.core_steady_state(parameters),
              parameters["beta"] / (parameters["k_g"] * parameters["a"]), 0.0)
    rhs = chemistry.deterministic_rhs(0.0, steady, parameters)
    np.testing.assert_allclose(rhs[:3], (0.0, 0.0, 0.0), atol=1e-14)
    assert rhs[3] == pytest.approx(parameters["beta"])
    eigenvalues = np.linalg.eigvals(chemistry.core_jacobian_at_steady(parameters))
    assert np.sign(max(value.real for value in eigenvalues)) == expected_sign
    assert chemistry.hopf_threshold(parameters["a"]) == 2.0


def test_generated_rules_parse_in_molfoundry():
    from mfobs.core.engine import mf
    from mfobs.core.gml import graph_gml, rule_gml

    for graph in [chemistry.species(name) for name in "AXYEBD"]:
        assert mf.Graph.fromGMLString(graph_gml(graph), add=False).isMolecule
    for rule in chemistry.rules().values():
        parsed = mf.Rule.fromGMLString(rule_gml(rule), add=False)
        assert parsed is not None


def test_public_molfoundry_application_reaches_declared_products():
    from mfobs.core.engine import mf
    from mfobs.core.gml import component_graphs, graph_gml, rule_gml

    cases = list(chemistry.core_rules().items())
    growth_host = disjoint_union([rop.chain(()), rop.monomer("B"), chemistry.species("X")])
    growth_target = disjoint_union([rop.chain(("B",)), chemistry.species("X")])
    cases.extend([
        (chemistry.MONOMER_INLET, chemistry.monomer_inlet_rule()),
        (chemistry.GROWTH, chemistry.growth_rule()),
    ])
    for name, rule in cases:
        host = growth_host if name == chemistry.GROWTH else rule_graph(rule, "left")
        expected = growth_target if name == chemistry.GROWTH else rule_graph(rule, "right")
        database, inputs = {}, []
        for part in component_graphs(host):
            parsed = mf.Graph.fromGMLString(graph_gml(part), add=False)
            parsed = database.setdefault(parsed.canonKey, parsed)
            inputs.append(parsed)
        engine_rule = mf.Rule.fromGMLString(rule_gml(rule), add=False)
        dg = mf.DG(graphDatabase=list(database.values()))
        with dg.build() as builder:
            outcomes = builder.apply(inputs, engine_rule, mode=mf.ApplyMode.UseAny, verbosity=0)
        assert len(outcomes) == 1, name
        actual_keys = sorted(vertex.graph.canonKey for vertex in outcomes[0].targets)
        expected_keys = sorted(
            mf.Graph.fromGMLString(graph_gml(part), add=False).canonKey
            for part in component_graphs(expected)
        )
        assert actual_keys == expected_keys, name


def test_exact_factorial_generator_finds_the_first_new_observable():
    assert moments.channel_derivative((1, 0), "source") == {
        (0, 0): moments.A_FEED * moments.VOLUME
    }
    assert moments.channel_derivative((1, 0), "autocatalysis") == {
        (2, 1): moments.VOLUME ** -2
    }
    assert moments.channel_derivative((1, 0), "conversion") == {
        (1, 0): -moments.B_FEED
    }
    assert moments.channel_derivative((1, 0), "sink") == {(1, 0): -1}


@pytest.mark.parametrize("moment", [(1, 0), (2, 1), (1, 1), (2, 2), (3, 1)])
@pytest.mark.parametrize("channel", sorted(moments.CHANNELS))
def test_factorial_conversion_matches_direct_generator_values(moment, channel):
    import sympy as sp

    count_x, count_y, volume = 7, 5, 11
    parameters = {
        moments.NX: count_x,
        moments.NY: count_y,
        moments.VOLUME: volume,
        moments.A_FEED: sp.Rational(6, 5),
        moments.B_FEED: sp.Rational(11, 5),
    }
    record = moments.CHANNELS[channel]
    dx, dy = record["delta"]
    observable = moments.factorial_expression(moment)
    direct = record["propensity"] * (
        observable.subs(
            {moments.NX: moments.NX + dx, moments.NY: moments.NY + dy},
            simultaneous=True,
        )
        - observable
    )
    expanded = sum(
        coefficient * moments.factorial_expression(target)
        for target, coefficient in moments.channel_derivative(moment, channel).items()
    )
    assert sp.cancel((direct - expanded).subs(parameters)) == 0


def test_budgeted_worklist_discovers_graph_moments_and_reports_unknown():
    from pathlib import Path

    from mfobs.core.worklist import run
    from mfobs.inputs import load

    root = Path(__file__).resolve().parents[1]
    inputs = load(root / "examples" / "brusselator-growth" / "discovery-model.json")
    result = run(inputs, max_expansions=4)
    assert result["closed"] is False
    assert result["status"] == "UNKNOWN"
    assert result["basis"][:4] == ["X", "Y", "ONE", "X2Y"]
    assert result["depths"]["X2Y"] == 1
    assert result["equations_by_reaction"]["X"]["autocatalysis"] == {
        "X2Y": "1/(1152*V**2)"
    }
    assert result["normalisation"]["X"]["scale"] == "24"
    assert result["normalisation"]["Y"]["scale"] == "48"
    assert result["normalisation"]["X2Y"]["scale"] == "27648"
    assert "X4Y" in result["basis"]
    assert result["moment_coordinates"]["X4Y"]["factorial_moment"] == "(N_X)_4 N_Y"
    assert result["raw_compositions"] == 0
    assert result["population_stoichiometry_from_gml"] == {
        "source": ({"A": 1}, {"X": 1}),
        "autocatalysis": ({"X": 2, "Y": 1}, {"X": 3}),
        "conversion": ({"B": 1, "X": 1}, {"D": 1, "Y": 1}),
        "sink": ({"X": 1}, {"E": 1}),
    }




def test_event_rate_scaling_keeps_all_six_channels():
    parameters = dict(a=2., b=3., k1=5., k2=7., k3=11., k4=13., beta=17., k_g=19.)
    assert chemistry.event_rates(parameters, volume=10) == {
        "k1*a*V": 100., "k2/V^2": .07, "k3*b": 33.,
        "k4": 13., "beta": 17., "k_g": 1.9,
    }
    with pytest.raises(ValueError, match="volume"):
        chemistry.event_rates(parameters, volume=0)




