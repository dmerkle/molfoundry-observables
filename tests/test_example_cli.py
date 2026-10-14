"""The artifact entry points use declared inputs and the production engines."""
from __future__ import annotations

import json
from pathlib import Path

import molfoundry as mf

from mfobs.query_input import from_file
from examples.tcs.inspect import export_gml, load_article_table, render_computed_table
from examples.tcs.replay import load_rules, run as run_tcs
from examples.tcs.table1 import compute_article_table
from examples.aba import run as aba_entry
from mfobs.example_cli import EXAMPLE_RUNS, _default_query_record


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def test_tcs_rules_are_loaded_from_visible_gml():
    rules = load_rules()
    assert set(rules) == {"aldol", "query_unconditioned", "query_conditioned"}
    for text in rules.values():
        assert mf.Rule.fromGMLString(text, add=False)
    assert "constrainNoEdge [ source 0 target 5 ]" in rules["aldol"]
    assert "constrainNoEdge" not in rules["query_unconditioned"]
    assert "constrainNoEdge [ source 0 target 5 ]" in rules["query_conditioned"]


def test_tcs_replay_distinguishes_the_transported_condition():
    result = run_tcs()
    assert result["status"] == "PASS"
    assert result["variants"]["conditioned"]["selected_net_commutator_coefficients"] == {
        "before": -1, "after": 1,
    }
    assert result["host_witness"]["input_match_counts"]["conditioned"] == {
        "before": 0, "after": 1,
    }
    assert result["host_witness"]["representation"].endswith("not a rewriting rule.")
    assert mf.Graph.fromGMLString(result["host_witness"]["gml"], add=False)


def test_tcs_article_table_is_computed_and_every_row_has_gml(tmp_path):
    oracle = load_article_table()
    assert oracle["status"] == "MANUAL_TRANSCRIPTION"
    table = compute_article_table()
    assert table["status"] == "COMPUTED_DIFFERS_FROM_PUBLISHED_TABLE"
    assert table["totals"] == {
        "computed": {"r_plus_after_r_id": 13, "r_id_after_r_plus": 18},
        "edge_reflecting_automatic": {
            "r_plus_after_r_id": 9, "r_id_after_r_plus": 18,
        },
        "published": {"r_plus_after_r_id": 12, "r_id_after_r_plus": 18},
    }
    assert len(table["rows"]) == 22
    assert table["semantics"]["phase"] == (
        "plain-rule composition before condition completion"
    )
    assert "Computed aldol composition table" in render_computed_table()

    result = run_tcs()
    summary = export_gml(result, tmp_path / "tcs-gml")
    assert summary["computed_table"]["rows"] == 22
    assert (tmp_path / "tcs-gml/computed-table/README.md").is_file()
    computed = json.loads((tmp_path / "tcs-gml/computed-table/table.json").read_text(encoding="utf-8"))
    for row in computed["rows"]:
        assert mf.Rule.fromGMLString(
            (tmp_path / "tcs-gml/computed-table" / row["gml"]).read_text(encoding="utf-8"), add=False
        )
    assert summary["selected"] == {
        "inputs": 3,
        "selected_composites": 4,
        "selected_commutator_terms": 2,
        "host_witness": "host-witness.gml",
    }
    manifest = json.loads((tmp_path / "tcs-gml/selected-commutators.json").read_text(encoding="utf-8"))
    assert manifest["terms"]["unconditioned"] == []
    assert sorted(term["coefficient"] for term in manifest["terms"]["conditioned"]) == [-1, 1]
    for term in manifest["terms"]["conditioned"]:
        assert mf.Rule.fromGMLString((tmp_path / "tcs-gml" / term["gml"]).read_text(encoding="utf-8"), add=False)


def test_ring_and_rop_inputs_are_declarative_and_present():
    ring = json.loads((ROOT / "examples/ring8/model.json").read_text(encoding="utf-8"))
    assert ring["initial_observables"] == ["P0"]
    assert all((ROOT / "examples/ring8" / rule["file"]).is_file()
               for rule in ring["rules"])
    query = json.loads((EXAMPLES / "aba/query.json").read_text(encoding="utf-8"))
    assert query["motif"] == "ABA"
    assert query["alphabet"] == ["A", "B"]
    assert query["gate"] is False


def test_every_paper_model_has_a_short_artifact_entry():
    assert set(EXAMPLE_RUNS) == {
        "tcs", "brusselator", "chlorination-guarded", "chlorination-ablated", "ring8",
        "aba", "rop", "remy", "growing-ring",
    }
    for name in ("brusselator", "chlorination", "growing_ring", "remy", "ring8", "aba", "molecular_rop"):
        assert (EXAMPLES / name / "run.py").is_file()
        assert (EXAMPLES / name / "README.md").is_file()


def test_common_query_files_recover_the_declared_queries():
    assert _default_query_record("chlorination-guarded") is None
    assert _default_query_record("chlorination-ablated") is None
    rop_query = from_file(EXAMPLES / "aba/query-gml.json")
    assert aba_entry.terminal.query_motifs(rop_query) == ("ABA",)
    remy_query = from_file(EXAMPLES / "remy/query.json")
    assert remy_query.names == ("P3",)
    from mfobs.models.remy import inputs as remy_trees
    assert remy_trees.query_order(remy_query) == 3
    growing_ring_query = from_file(EXAMPLES / "growing_ring/query.json")
    assert growing_ring_query.names == ("Q0", "Q1", "Q2")
    assert _default_query_record("growing-ring") == growing_ring_query


def test_aba_profile_uses_the_same_question_as_the_artifact_input():
    query = json.loads(aba_entry.QUERY.read_text(encoding="utf-8"))
    profile = json.loads((ROOT / "scripts/reproduction/paper_evidence_manifest.json").read_text(encoding="utf-8"))
    stage = profile["stages"]["aba"]
    assert query == {
        "motif": stage["motif"], "alphabet": stage["alphabet"],
        "gate": stage["gated"], "linear_limit": stage["linear_limit"],
        "support_pattern_limit": stage["support_pattern_limit"],
        "support_expansion_limit": stage["support_expansion_limit"],
    }
