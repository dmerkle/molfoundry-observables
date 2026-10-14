import mfobs
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.core.composition_jump_model import compare_polynomial_jump_models
"""All-size construction and finite graph regressions for the ROP family."""
from mfobs.checks.terminal_paths import validate_result
from mfobs.applications import ABA

import mfobs.checks.molecular_rop as model_molecular_rop
import mfobs.checks.rop_projection as model_rop_projection
import json

import pytest

from mfobs.core.conditions import engine_no_edges, matches
from mfobs.core.graph import apply_match, edge_rows, node_labels
from mfobs.core.rules import added_edges
import molfoundry as mf
from mfobs.models.molecular_rop import inputs as model


def test_four_variants_share_one_normalized_analysis_contract():
    for name, selected in model.VARIANTS.items():
        contract = model_rop_projection.analysis_contract(name)
        grammar = model_rop_projection.grammar(name)
        assert contract.digest == grammar.input_definition_sha256
        assert contract.state_language.infinite
        assert all(rule.semantics == "conditional-DPO" for rule in contract.rules)
        assert contract.linear_seeds == tuple(
            observable.name
            for observable in grammar.primary_observables()
            if observable.role == "selected"
        )
        assert contract.pjc_targets == contract.linear_seeds
        assert all(
            observable.kind in {"symbol-count", "factor-count"}
            for observable in grammar.primary_observables()
            if observable.name in contract.linear_seeds
        )
        assert selected.purpose


def test_support_coordinates_are_not_silently_used_as_linear_seeds():
    backend1 = ABA(('B',), alphabet=('B',), gated=False)
    baseline = mfobs.analyze(backend1.specification(), backend=backend1).record()
    assert baseline["contract"]["analysis_branches"][
        "linear_observable_discovery"
    ]["seeds"] == ["N_B"]
    assert baseline["linear"]["origins"]["M_B"] == {
        "origin": "computed-channel-residual",
        "discovered_from": "N_B",
        "rate_channel": "k_B",
    }

    backend2 = ABA(('B',), alphabet=('B',), gated=True)
    noedge = mfobs.analyze(backend2.specification(), backend=backend2).record()
    assert noedge["contract"]["analysis_branches"][
        "linear_observable_discovery"
    ]["seeds"] == ["N_B"]
    assert noedge["linear"]["basis"] == ["N_B", "M_B*X", "X", "M_B", "1"]
    assert noedge["linear"]["origins"]["M_B*X"]["discovered_from"] == "N_B"
    roles = {
        observable["name"]: observable["role"]
        for observable in noedge["contract"]["observables"]
    }
    assert roles == {"N_B": "selected", "M_B": "auxiliary", "X": "auxiliary"}


def test_linear_discovery_and_pjc_are_distinct_analysis_branches():
    expected = {
        "baseline": ("CLOSED", 3, 1),
        "noedge": ("CLOSED", 5, 2),
        "copolymer": ("UNKNOWN_DIMENSION_LIMIT", 24, 2),
        "combined": ("UNKNOWN_DIMENSION_LIMIT", 24, 3),
    }
    for name, (linear_status, dimension, order) in expected.items():
        selected = model.configuration(name)
        backend4 = ABA(selected.alphabet + tuple(''.join(w) for w in selected.motifs), alphabet=selected.alphabet, gated=selected.gated)
        result = mfobs.analyze(backend4.specification(), backend=backend4, linear_limit=24).record()
        assert result["linear"]["status"] == linear_status
        assert len(result["linear"]["basis"]) == dimension
        assert result["pjc"]["status"] == "CLOSED"
        assert result["pjc"]["operator"]["maximum_differential_order"] == order
        assert result["linear"]["input_definition_sha256"] == result["contract"]["definition_sha256"]
        assert result["pjc"]["input_definition_sha256"] == result["contract"]["definition_sha256"]


def test_all_presets_use_query_driven_support_and_uniform_event_extraction():
    expected_classes = {"baseline": 2, "noedge": 4, "copolymer": 10, "combined": 12}
    for name, count in expected_classes.items():
        selected = model.configuration(name)
        backend5 = ABA(selected.alphabet + tuple(''.join(w) for w in selected.motifs), alphabet=selected.alphabet, gated=selected.gated)
        result = mfobs.analyze(backend5.specification(), backend=backend5).record()
        evidence = result["composition"]
        assert evidence["status"] == validate_result(result)["status"] == "PASS"
        assert evidence["jump_projection"] == "DPO"
        assert len(evidence["event_classes"]) == count
        assert evidence["projection"]["structural_class_count"] == count
        assert set(evidence["observable_coordinates"]) == set(result["support_discovery"]["patterns"])
        assert "independent_word_oracle" not in evidence
        provenance = evidence["jump_model_provenance"]
        assert provenance["all_jumps_composition_derived"] is True
        assert provenance["supplied_delta_table"] is False


def test_structural_projection_preserves_matches_targets_and_observables():
    for name in model.VARIANTS:
        result = model_molecular_rop.composition_projection_validation(
            name, maximum_word_length=2
        )
        assert result["status"] == "PASS", result["failures"][:1]
        assert result["observable_match_checks"]
        assert result["structural_class_match_checks"]
        assert result["structural_target_checks"]
        assert result["molecular_to_structural_multiplicity_checks"]
        assert not result["finite_enumeration_used_by_analysis"]


def test_suffix_automaton_is_generated_from_selected_motif_depth():
    grammar = model_rop_projection.grammar("copolymer")
    assert grammar.suffix_depth == 2
    assert ["".join(value) or "empty" for value in grammar.suffix_classes] == [
        "empty", "A", "B", "AA", "AB", "BA", "BB"
    ]
    auxiliary = {value.name for value in grammar.auxiliary_observables()}
    assert auxiliary == {"S_A", "S_B", "S_AA", "S_AB", "S_BA", "S_BB"}


@pytest.mark.parametrize("gated", [False, True])
def test_custom_motif_query_generates_its_own_suffix_partition_and_gate(gated):
    custom = model.configured_variant(motifs=("ABBA",), gated=gated)
    language = model_rop_projection.grammar(custom)
    graph_model = TerminalWordGraphModel(language)

    assert language.suffix_depth == 3
    assert language.motifs == (("A", "B", "B", "A"),)
    assert len(graph_model.reaction_classes()) == 32 + (2 if gated else 0)
    assert graph_model.event_partition_certificate()["suffix_depth"] == 3
    assert {
        "append-A-from-ABB", "append-B-from-ABB"
    } <= {event.identifier for event in graph_model.reaction_classes()}
    if gated:
        rules = model.rules(custom)
        open_state = model.state_graph(custom, model.WordState((), (1, 1), True))
        closed_state = model.state_graph(custom, model.WordState((), (1, 1), False))
        assert {"donor_bind", "donor_unbind"} <= set(rules)
        assert len(matches(rules["insert_A"], open_state)) == 1
        assert len(matches(rules["insert_A"], closed_state)) == 0
    else:
        assert all(
            set(engine_no_edges(rule)) <= {(u, v) for u, v, _ in added_edges(rule)}
            for rule in model.rules(custom).values()
        )


def test_custom_motif_query_is_accepted_by_the_composition_compiler():
    # A non-preset query exercises the composition route itself; the ABBA
    # regression above additionally exercises the deeper generated refinement.
    custom = model.configured_variant(motifs=("ABB",), gated=False)
    compilation = model_rop_projection.composition_compilation(custom)
    oracle = model_rop_projection.grammar(custom).polynomial_jump_model()

    assert len(compilation.model.events) == 16
    assert "C_ABB" in compilation.model.ring.variables
    assert compare_polynomial_jump_models(compilation.model, oracle)["status"] == "PASS"


def test_generated_refinement_check_covers_exact_sources_and_rejects_a_gap():
    query = model.configured_variant(
        motifs=("AA", "AB", "BA", "BB", "ABBA"), gated=True,
    )
    projection = TerminalWordGraphModel(model_rop_projection.grammar(query))
    classes = projection.reaction_classes()
    passed = projection.refinement_consistency_check(classes)
    incomplete = [
        reaction for reaction in classes
        if reaction.identifier != "append-A-from-ABB"
    ]
    failed = projection.refinement_consistency_check(incomplete)

    assert passed["status"] == "PASS"
    assert passed["maximum_motif_length"] == 4
    assert len(passed["expected_source_rules"]["insert_A"]) == 15
    assert failed["status"] == "FAIL"
    assert failed["failures"][0]["source_rule"] == "insert_A"


def test_preset_configuration_remains_reproducible_and_budget_is_hashed():
    first = model.effective_configuration("copolymer", linear_limit=24)
    second = model.effective_configuration("copolymer", linear_limit=25)
    short_query = model.configured_variant(motifs=("ABB",), gated=False)
    long_query = model.configured_variant(motifs=("ABBA",), gated=False)
    backend3 = ABA(('A', 'B', 'AA', 'AB', 'BA', 'BB', 'ABA'), alphabet=('A', 'B'), gated=False)
    result = mfobs.analyze(backend3.specification(), backend=backend3, linear_limit=24).record()

    assert first["model"] == model.variant("copolymer").record()
    assert first["sha256"] != second["sha256"]
    assert short_query.identifier != long_query.identifier
    assert (
        model_rop_projection.analysis_contract(short_query).digest
        != model_rop_projection.analysis_contract(long_query).digest
    )
    assert result["configuration"]["model"]["motifs"] == ["A", "B", "AA", "AB", "BA", "BB", "ABA"]
    assert result["configuration"]["model"]["alphabet"] == first["model"]["alphabet"]
    assert result["configuration"]["model"]["gate"] == first["model"]["gate"]
    assert result["configuration"]["linear_dimension_limit"] == first["linear_dimension_limit"]
    assert result["configuration"]["representation"] == "directed paths"
    assert len(result["composition"]["event_classes"]) == 10


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"motifs": ("AC",), "gated": False}, "outside the selected alphabet"),
        ({"motifs": ("AB", "AB"), "gated": False}, "must be unique"),
        ({"motifs": ("A",), "gated": False}, "length at least two"),
    ],
)
def test_custom_configuration_rejects_invalid_motif_queries(arguments, message):
    with pytest.raises(ValueError, match=message):
        model.configured_variant(**arguments)


def test_configuration_rejects_nonpositive_linear_budget():
    with pytest.raises(ValueError, match="positive integer"):
        model.effective_configuration("copolymer", linear_limit=0)


def test_too_small_linear_budget_fails_before_composition(monkeypatch):
    query = model.configured_variant(
        motifs=("AA", "AB", "BA", "BB", "ABBA"), gated=False,
    )
    compiled = False

    def must_not_compile(_selection):
        nonlocal compiled
        compiled = True
        raise AssertionError("composition must not run for an impossible budget")

    from mfobs.models.aba import counting as terminal_paths
    monkeypatch.setattr(terminal_paths, "derive", must_not_compile)
    with pytest.raises(ValueError, match="smaller than the 5 selected"):
        backend6 = ABA(tuple((''.join(word) for word in query.motifs)) or query.alphabet, alphabet=query.alphabet, gated=query.gated)
        mfobs.analyze(backend6.specification(), backend=backend6, linear_limit=4).record()
    assert not compiled


def test_cli_exposes_a_nonreport_custom_query_entry_point():
    from mfobs.cli import parser

    arguments = parser().parse_args([
        "rop-analysis", "--motif", "ABBA", "--gate", "--linear-limit", "17",
        "--out", "/tmp/rop-abba.json",
    ])
    assert arguments.command == "rop-analysis"
    assert arguments.motif == ["ABBA"]
    assert arguments.gate is True
    assert arguments.linear_limit == 17


def test_cli_rejects_an_explicitly_empty_alphabet(capsys, tmp_path):
    from mfobs.cli import main

    status = main([
        "rop-analysis", "--motif", "AB", "--alphabet", "",
        "--out", str(tmp_path / "unused.json"),
    ])
    assert status == 1
    assert "provided but names no chemical symbols" in capsys.readouterr().err


def test_cli_reports_bounded_statuses_events_and_diagnostic_timings(capsys, tmp_path):
    from mfobs.cli import main

    out = tmp_path / "custom.json"
    status = main([
        "rop-analysis", "--motif", "ABB", "--gate", "--linear-limit", "24",
        "--diagnostic-timings", "--out", str(out),
    ])
    assert status == 0
    message = capsys.readouterr().out
    assert "events=10" in message
    assert "linear=UNKNOWN_DIMENSION_LIMIT" in message
    assert "pJC=CLOSED" in message
    assert "oracle=" not in message
    assert "Scope: exact analysis is relative to the declared ROP state language" in message
    result = json.loads(out.read_text(encoding="utf-8"))
    assert result["diagnostics"]["seconds"]
    assert all(value >= 0 for value in result["diagnostics"]["seconds"].values())


def test_noedge_gate_is_active_and_not_an_inert_annotation():
    selected = model.variant("combined")
    grammar = model_rop_projection.grammar("combined")
    open_state = model.state_graph(
        "combined", model.WordState((), (1, 1), True)
    )
    closed_state = model.state_graph(
        "combined", model.WordState((), (1, 1), False)
    )
    for symbol in selected.alphabet:
        rule = model.rules("combined")[f"insert_{symbol}"]
        assert engine_no_edges(rule)
        mappings = matches(rule, open_state)
        assert len(mappings) == 1
        assert not matches(rule, closed_state)
        product = apply_match(rule, open_state, mappings[0])
        expected = model.state_graph("combined", model.WordState(
            (symbol,), tuple(int(other != symbol) for other in selected.alphabet), True
        ))
        assert product is not None and product.canonKey == expected.canonKey
    assert grammar.gated


def test_lactone_insertion_is_atom_preserving_conditional_dpo():
    for symbol in ("A", "B"):
        left, right = model.propagation_sides(symbol)
        assert node_labels(left) == node_labels(right)
        assert sorted(node_labels(left).values()) == sorted(node_labels(right).values())
        assert len(set(edge_rows(left)) - set(edge_rows(right))) == 2
        assert len(set(edge_rows(right)) - set(edge_rows(left))) == 2


def test_literal_graph_observables_and_rule_matches_agree_with_word_semantics():
    for name in model.VARIANTS:
        result = model_molecular_rop.graph_semantics_validation(name, maximum_word_length=2)
        assert result["status"] == "PASS", result["failures"][:1]
        assert not result["finite_enumeration_used_by_analysis"]


def test_generated_pjc_matches_direct_semantic_events():
    for name in model.VARIANTS:
        result = model_molecular_rop.polynomial_operator_validation(name, maximum_word_length=3)
        assert result["status"] == "PASS", result["failures"][:1]
        assert not result["finite_enumeration_used_by_analysis"]


def test_gml_exports_conditions_and_every_observable():
    exported = model.exported_gml("combined")
    assert set(exported["rules"]) == {
        "feed_A", "feed_B", "insert_A", "insert_B", "donor_bind", "donor_unbind"
    }
    native = model.rules("combined")
    for name in ("insert_A", "donor_bind"):
        restored = mf.Rule.fromGMLString(exported["rules"][name], add=False)
        assert engine_no_edges(restored) == engine_no_edges(native[name])
        assert restored.canonKey == native[name].canonKey
    # MolFoundry need not serialize created-edge admissibility as an extra guard.
    assert (0, 4) in engine_no_edges(native["donor_bind"])
    assert {"C_AA", "C_AB", "C_BA", "C_BB", "C_ABA"} <= set(exported["observables"])
