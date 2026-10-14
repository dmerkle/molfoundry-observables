"""Standard CLI/API integration of discovered support with existing consumers."""
from mfobs.checks.terminal_paths import validate_result
import mfobs
from mfobs.applications import ABA

import mfobs.models.aba.events as discovery_events
import mfobs.models.aba.inputs as discovery_inputs
import mfobs.checks.rop_projection as rop_rop_projection
from dataclasses import replace
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys

import pytest

from mfobs.core.pattern_support import UnsupportedReduction, polynomial_record
from mfobs.core.polynomial_jump_closure import generator_action
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.aba import counting as discovery


def custom(word="ABA", gated=False):
    return rop.configured_variant(motifs=(word,), gated=gated)


@pytest.mark.parametrize("gated", [False, True])
def test_only_aba_is_selected_and_existing_consumers_receive_discovered_events(gated):
    backend1 = ABA(tuple((''.join(word) for word in custom(gated=gated).motifs)) or custom(gated=gated).alphabet, alphabet=custom(gated=gated).alphabet, gated=custom(gated=gated).gated)
    result = mfobs.analyze(backend1.specification(), backend=backend1, linear_limit=24).record()
    selected = result["contract"]["analysis_branches"]["linear_observable_discovery"]["seeds"]
    assert selected == ["C_ABA"]
    assert result["support_discovery"]["selected_queries"] == selected
    assert set(result["pjc"]["observable_coordinates"]) == {
        "C_ABA", "M_A", "M_B", "S_AB", "S_A", *(["X"] if gated else []),
    }
    assert result["support_discovery"]["patterns"]["S_AB"]["origin"]["pattern"] == "C_ABA"
    assert result["support_discovery"]["patterns"]["S_A"]["origin"]["pattern"] == "S_AB"
    assert result["linear"]["status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert len(result["linear"]["basis"]) == 24
    assert result["pjc"]["status"] == "CLOSED"
    assert result["pjc"]["event_class_count"] == 8 + 2 * gated
    assert result["pjc"]["operator"]["maximum_differential_order"] == 2 + gated
    assert validate_result(result)["status"] == "PASS"
    for branch in ("linear", "pjc"):
        assert result[branch]["input_definition_sha256"] == result["contract"]["definition_sha256"]
    assert result["pjc"]["jump_model_provenance"]["supplied_delta_table"] is False


def test_second_motif_is_discovered_without_building_a_full_suffix_grammar(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("discovery must not enumerate a full suffix basis")
    monkeypatch.setattr(rop, "TerminalGrowthGrammar", forbidden)
    monkeypatch.setattr(rop_rop_projection, "composition_compilation", forbidden)
    backend2 = ABA(tuple((''.join(word) for word in custom('ABBA', gated=True).motifs)) or custom('ABBA', gated=True).alphabet, alphabet=custom('ABBA', gated=True).alphabet, gated=custom('ABBA', gated=True).gated)
    result = mfobs.analyze(backend2.specification(), backend=backend2, linear_limit=16).record()
    patterns = result["support_discovery"]["patterns"]
    assert set(patterns) == {"C_ABBA", "M_A", "M_B", "X", "S_ABB", "S_AB", "S_A"}
    assert result["support_discovery"]["selected_queries"] == ["C_ABBA"]
    assert validate_result(result)["status"] == "PASS"
    assert result["pjc"]["operator"]["maximum_differential_order"] == 3


def test_adapter_preserves_each_derived_first_generator_image_and_suffix_identities():
    support, sources = discovery.derive(("ABA",), gated=True)
    model = discovery_events.polynomial_jump_model(support, sources, name="test", input_definition_sha256="test")
    ring = model.ring
    assert ring.monomial({"S_A": 1, "S_AB": 1}).is_zero
    for name in ring.variables:
        for source in sources:
            symbolic = polynomial_record(support.images[name][source.identifier])
            expected = ring.add(*[
                ring.scale(ring.monomial(t["powers"]), Fraction(t["coefficient"]))
                for t in symbolic["terms"]
            ])
            assert generator_action(model, ring.coordinate(name), channel=source.rate) == expected
    # A different query supplies nested, compatible suffix indicators too.
    other, rules = discovery.derive(("AAAA",))
    longer = discovery_events.polynomial_jump_model(other, rules, name="nested", input_definition_sha256="test")
    assert longer.ring.monomial({"S_A": 1, "S_AAA": 1}) == longer.ring.coordinate("S_AAA")


@pytest.mark.parametrize("budget,status", [
    ({"support_pattern_limit": 3}, "UNKNOWN_PATTERN_LIMIT"),
    ({"support_expansion_limit": 0}, "UNKNOWN_EXPANSION_LIMIT"),
])
def test_support_budgets_preserve_pending_obligations_and_do_not_export_an_operator(budget, status):
    backend3 = ABA(tuple((''.join(word) for word in custom().motifs)) or custom().alphabet, alphabet=custom().alphabet, gated=custom().gated)
    result = mfobs.analyze(backend3.specification(), backend=backend3, **budget).record()
    assert result["status"] == status
    assert result["support_discovery"]["pending"]
    assert result["support_discovery"]["events"]["status"] == "UNKNOWN_SUPPORT"
    assert validate_result(result)["status"] == "NOT_EXECUTED"
    assert result["linear"]["status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert result["pjc"]["status"] == "UNKNOWN_PREREQUISITE"
    assert "operator" not in result["pjc"]


def test_unsupported_reduction_is_not_silently_replaced_by_the_suffix_route(monkeypatch):
    def unsupported(*args):
        raise UnsupportedReduction("test unsupported graph condition")
    monkeypatch.setattr(discovery.PathCounter, "image", unsupported)
    backend4 = ABA(tuple((''.join(word) for word in custom().motifs)) or custom().alphabet, alphabet=custom().alphabet, gated=custom().gated)
    result = mfobs.analyze(backend4.specification(), backend=backend4).record()
    assert result["status"] == "UNKNOWN_UNSUPPORTED_REDUCTION"
    assert result["support_discovery"]["reason"] == "test unsupported graph condition"
    assert result["composition"]["event_classes"] == []


def test_event_uniformity_failure_stays_separate_from_successful_support(monkeypatch):
    original = discovery_events.compile_uniform_events
    def failed(support, sources):
        if support.status == "POLYNOMIAL_SUPPORT_CLOSED":
            return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "events": []}
        return original(support, sources)
    monkeypatch.setattr(discovery_events, "compile_uniform_events", failed)
    backend5 = ABA(tuple((''.join(word) for word in custom().motifs)) or custom().alphabet, alphabet=custom().alphabet, gated=custom().gated)
    result = mfobs.analyze(backend5.specification(), backend=backend5).record()
    assert result["support_discovery"]["status"] == "POLYNOMIAL_SUPPORT_CLOSED"
    assert result["status"] == "UNKNOWN_UNIFORM_MATCH_OBLIGATION"
    assert result["linear"]["status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert result["pjc"]["status"] == "UNKNOWN_PREREQUISITE"


@pytest.mark.parametrize("preset", ["baseline", "noedge"])
def test_named_presets_discover_support_and_preserve_exact_linear_dynamics(preset):
    backend6 = ABA(("B",), alphabet=("B",), gated=preset == "noedge")
    result = mfobs.analyze(backend6.specification(), backend=backend6, linear_limit=24).record()
    backend9 = ABA(("B",), alphabet=("B",), gated=preset == "noedge")
    assert result == mfobs.analyze(backend9.specification(), backend=backend9, linear_limit=24).record()
    assert result["support_discovery"]["selected_queries"] == ["N_B"]
    assert result["linear"]["status"] == result["pjc"]["status"] == "CLOSED"
    # Independent full suffix refinement; the discovered model need not use its basis.
    reference = rop_rop_projection.composition_compilation(preset).model
    ring = reference.ring
    def decode(record):
        return ring.add(*[ring.scale(ring.monomial(t['monomial']), Fraction(t['coefficient']))
                          for t in record['terms']])
    linear = result['linear']
    expressions = {n: decode(v) for n,v in linear['basis_expressions'].items()}
    for step in linear['steps']:
        actual = generator_action(reference, expressions[step['source']], channel=step['channel'])
        reconstructed = ring.add(*[ring.scale(expressions[n], Fraction(c))
                                   for n,c in step['coordinates'].items()])
        assert actual == decode(step['generator_image']) == reconstructed
    assert expressions['N_B'] == ring.coordinate('N_B')


def test_applications_do_not_accept_an_implicit_representation_switch():
    with pytest.raises(TypeError, match="support_method"):
        ABA.motif("ABA", support_method="molecular")


def test_standard_command_writes_discovery_json_and_reports_unknown_without_pass(tmp_path):
    root = Path(__file__).resolve().parents[1]
    out = tmp_path / "analysis.json"
    base = [sys.executable, "run.py", "rop-analysis", "--motif", "ABA", "--out", str(out)]
    completed = subprocess.run(base, cwd=root, capture_output=True, text=True, check=True)
    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["support_discovery"]["selected_queries"] == ["C_ABA"]
    assert "events=8" in completed.stdout
    limited = subprocess.run(base + ["--support-expansion-limit", "0"], cwd=root,
                             capture_output=True, text=True, check=True)
    assert "UNKNOWN_EXPANSION_LIMIT:" in limited.stdout
    assert "PASS:" not in limited.stdout
    assert json.loads(out.read_text(encoding="utf-8"))["pjc"]["status"] == "UNKNOWN_PREREQUISITE"


def test_invalid_budgets_are_rejected_before_composition(monkeypatch):
    monkeypatch.setattr(discovery, "derive", lambda *a, **k: pytest.fail("composition ran"))
    for kwargs in ({"support_expansion_limit": -1}, {"support_pattern_limit": 0},
                   {"support_pattern_limit": True}, {"linear_limit": 0}):
        with pytest.raises(ValueError):
            backend11 = ABA(tuple((''.join(word) for word in custom().motifs)) or custom().alphabet, alphabet=custom().alphabet, gated=custom().gated)
            mfobs.analyze(backend11.specification(), backend=backend11, **kwargs).record()


def test_source_metadata_and_unsupported_graph_shapes_cannot_be_coerced():
    support, sources = discovery.derive(("ABA",))
    changed = (replace(sources[0], kind="append"), *sources[1:])
    assert discovery_events.compile_uniform_events(support, changed)["status"] == "UNKNOWN_UNIFORM_MATCH_OBLIGATION"
    bad = [
        discovery_inputs.graph_rule({0: "Feed:A"}, [], adjacency=[(0, ">", 0)]),
        discovery_inputs.graph_rule({0: "Position", 1: "Symbol:A", 2: "Ghost"},
                             [(0, 1, "kind"), (0, 2, "annotation")]),
        discovery_inputs.graph_rule({0: "Position", 1: "Symbol:A", 2: "Position", 3: "Symbol:B"},
                             [(0, 1, "kind"), (2, 3, "kind"), (1, 3, "annotation")]),
    ]
    for rule in bad:
        with pytest.raises(UnsupportedReduction):
            discovery.reduce_input(rule)


def test_optional_composition_record_preserves_analysis_and_signed_multiplicities():
    from collections import Counter
    import molfoundry as mf

    backend7 = ABA(tuple((''.join(word) for word in custom(gated=True).motifs)) or custom(gated=True).alphabet, alphabet=custom(gated=True).alphabet, gated=custom(gated=True).gated)
    baseline = mfobs.analyze(backend7.specification(), backend=backend7, linear_limit=24).record()
    backend8 = ABA(tuple((''.join(word) for word in custom(gated=True).motifs)) or custom(gated=True).alphabet, alphabet=custom(gated=True).alphabet, gated=custom(gated=True).gated)
    recorded = mfobs.analyze(backend8.specification(), backend=backend8, linear_limit=24, record_compositions=True).record()
    saw_cancellation = saw_rejection = False
    for step in recorded["support_discovery"]["trace"]:
        if "composition_orders" not in step:
            continue
        products = step.pop("raw_composites")
        classes = step.pop("conditional_classes")
        coefficients = Counter()
        for order, diagnostic in step["composition_orders"].items():
            rows = [p for p in products if p["order"] == order]
            assert len(rows) == diagnostic["public_raw"]
            assert sum(not p["rejection_reason"] for p in rows) == diagnostic["retained"]
        for row in products:
            rule = mf.Rule.fromGMLString(row["gml"], add=False)
            if row["rejection_reason"]:
                saw_rejection = True
            else:
                coefficients[rule.canonKey] += row["coefficient"]
        assert {
            mf.Rule.fromGMLString(row["gml"], add=False).canonKey: row["coefficient"]
            for row in classes
        } == dict(coefficients)
        saw_cancellation |= any(value == 0 for value in coefficients.values())
    assert saw_cancellation and saw_rejection
    assert recorded == baseline


@pytest.mark.parametrize("method", ["discovery", "molecular"])
def test_analysis_does_not_run_reference_checks(monkeypatch, method):
    from mfobs.checks import terminal_paths, molecular_rop
    def forbidden(*args, **kwargs):
        raise AssertionError("reference comparison called by analysis")
    monkeypatch.setattr(terminal_paths, "validate_words", forbidden)
    monkeypatch.setattr(molecular_rop, "validate_molecular_hosts", forbidden)
    monkeypatch.setattr(molecular_rop, "compare_typed_reference", forbidden)
    backend = ABA(('B',), alphabet=('B',), gated=False)
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    assert result["status"] == "ANALYZED"
    assert not set(result) & {"word_oracle_comparison", "molecular_graph_validation", "typed_reference_comparison"}


def test_external_word_check_detects_a_changed_exported_joint_jump():
    from copy import deepcopy
    backend = ABA.motif('ABA')
    result = mfobs.analyze(backend.specification(), backend=backend).record()
    assert validate_result(result, maximum_length=3)["status"] == "PASS"
    damaged = deepcopy(result)
    damaged["composition"]["event_classes"][0]["jump"]["C_ABA"] = 99
    with pytest.raises(AssertionError):
        validate_result(damaged, maximum_length=3)
