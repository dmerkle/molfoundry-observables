import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from scripts.reproduction import reproduce_paper as reproduction


def test_manifest_profiles_reference_declared_stages():
    manifest = reproduction.load_manifest()
    assert set(manifest["profiles"]) == {"smoke", "paper", "supplementary"}
    required = {"doctor", "aldol-conditions", "molecular-baseline", "molecular-limited-aba"}
    assert required < set(manifest["profiles"]["smoke"])
    assert required < set(manifest["profiles"]["paper"])
    assert {"context-nac-guarded", "context-nac-ablated"} <= set(manifest["profiles"]["paper"])
    for stage in manifest["stages"].values():
        for test in stage.get("tests", []):
            assert (reproduction.PROJECT / test.split("::", 1)[0]).is_file(), test
    assert len(set(manifest["profiles"]["paper"])) == len(manifest["profiles"]["paper"])


def test_manifest_rejects_foreign_schema(tmp_path):
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps({"schema_version": 2, "kind": "foreign", "stages": {}}), encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported"):
        reproduction.load_manifest(path)


def test_expected_unknown_and_partial_are_exact_structured_assertions():
    growing = {"linear_status": "UNKNOWN_EXPANSION_LIMIT",
               "pjc_status": "UNKNOWN_UNBOUNDED_COORDINATE_LIFT",
               "witness_status": "WITNESS"}
    assert reproduction.compare_expected(growing, growing)["status"] == "PASS"
    partial = {"status": "COMPARISON_INCOMPLETE", "minimum_seed_equalities": 7,
               "mapped_mismatches": 0, "finite_status": "PASS"}
    assert reproduction.compare_expected(partial, partial)["status"] == "PASS"
    partial["mapped_mismatches"] = 1
    assert reproduction.compare_expected(partial, {
        **partial, "mapped_mismatches": 0})["status"] == "FAIL"


def test_missing_expected_field_is_not_silently_accepted():
    checked = reproduction.compare_expected({"status": "PASS"},
                                            {"status": "PASS", "event_count": 4})
    assert checked["status"] == "FAIL"
    assert checked["fields"]["event_count"]["actual"] is None


def test_support_discovery_stage_is_the_single_aba_query():
    stage = reproduction.load_manifest()["stages"]["aba"]
    assert stage["worker"] == "aba"
    assert stage["motif"] == "ABA"
    assert stage["alphabet"] == ["A", "B"]
    assert stage["gated"] is False
    assert stage["expected"]["event_count"] == 8
    assert stage["expected"]["support_status"] == "POLYNOMIAL_SUPPORT_CLOSED"


def test_tcs_stage_keeps_local_replay_scope():
    stage = reproduction.load_manifest()["stages"]["tcs-aldol"]
    assert stage["worker"] == "tcs-aldol"
    assert stage["expected"]["conditioned_coefficients"] == [-1, 1]
    assert stage["expected"]["conditioned_host_counts"] == [0, 1]


def test_input_identity_requires_every_declared_file(tmp_path):
    present = tmp_path / "present.json"
    present.write_text("{}\n", encoding="utf-8")
    identity = reproduction.input_identity([present])
    assert list(identity["files"]) == [str(present.resolve())]
    assert identity["files"][str(present.resolve())] == reproduction.sha256(present)
    with pytest.raises(FileNotFoundError, match="required input evidence"):
        reproduction.input_identity([present, tmp_path / "missing.json"])


def test_rop_worker_rejects_checks_for_a_different_input(monkeypatch):
    from mfobs.checks import molecular_rop, rop_projection

    import mfobs
    monkeypatch.setattr(mfobs, "analyze", lambda *args, **kwargs: type("Result", (), {"record": lambda self: {}})())
    monkeypatch.setattr(rop_projection, "analysis_contract",
                        lambda *args: SimpleNamespace(digest="expected"))
    for name in ("graph_semantics_validation", "composition_projection_validation",
                 "polynomial_operator_validation"):
        monkeypatch.setattr(molecular_rop, name,
                            lambda *args: {"input_definition_sha256": "foreign"})
    with pytest.raises(ValueError, match="different ROP input"):
        reproduction.worker_rop("rop-baseline", {
            "selection": "baseline", "linear_limit": 24, "graph_bound": 1,
            "projection_bound": 1, "operator_bound": 1})


def test_failed_child_process_is_retained_as_process_error(tmp_path):
    record = reproduction.run_process(
        [sys.executable, "-c", "raise SystemExit(7)"], stage="simulated-failure",
        output=tmp_path, timeout=10,
    )
    assert record["execution_status"] == "PROCESS_ERROR"
    assert record["process_returncode"] == 7
    assert Path(record["stderr_log"]).is_file()


def test_verbose_child_is_drained_directly_to_separate_logs(tmp_path):
    record = reproduction.run_process(
        [sys.executable, "-c",
         "import sys;sys.stdout.write('x'*2000000);sys.stderr.write('y'*2000000)"],
        stage="verbose-child", output=tmp_path, timeout=10,
    )
    assert record["execution_status"] == "COMPLETE"
    assert Path(record["stdout_log"]).stat().st_size == 2_000_000
    assert Path(record["stderr_log"]).stat().st_size == 2_000_000


def test_failed_worker_cannot_reuse_stale_success_record(tmp_path, monkeypatch):
    manifest = reproduction.load_manifest()
    stage_dir = tmp_path / "stages" / "doctor"
    stage_dir.mkdir(parents=True)
    stale = {"schema_version": 1, "kind": "mfobs-stage-worker",
             "stage": "doctor", "execution_status": "COMPLETE",
             "reproduction_status": "EXPECTED", "worker_elapsed_seconds": 1,
             "expectation_check": {"status": "PASS"}, "observed": {"semantic_status": "PASS"},
             "inputs": {}}
    (stage_dir / "record.json").write_text(json.dumps(stale), encoding="utf-8")
    failing = tmp_path / "failing.py"
    failing.write_text("raise SystemExit(7)\n", encoding="utf-8")
    monkeypatch.setattr(reproduction, "__file__", str(failing))
    result = reproduction._run_stage(
        "doctor", manifest["stages"]["doctor"],
        SimpleNamespace(manifest=reproduction.MANIFEST_PATH), tmp_path,
    )
    assert result["process_returncode"] == 7
    assert result["execution_status"] == "PROCESS_ERROR"
    assert result["reproduction_status"] == "UNEXPECTED"
    assert result["outcome_class"] == "FAILED"


@pytest.mark.parametrize("payload", ["{", '{"schema_version": 1}'])
def test_malformed_or_incomplete_worker_record_is_invalid_output(tmp_path, monkeypatch, payload):
    manifest = reproduction.load_manifest()
    emitter = tmp_path / "emitter.py"
    emitter.write_text(
        "import pathlib, sys\npathlib.Path(sys.argv[-1]).write_text(" + repr(payload) + ")\n", encoding="utf-8"
    )
    monkeypatch.setattr(reproduction, "__file__", str(emitter))
    result = reproduction._run_stage(
        "doctor", manifest["stages"]["doctor"],
        SimpleNamespace(manifest=reproduction.MANIFEST_PATH), tmp_path,
    )
    assert result["execution_status"] == "INVALID_OUTPUT"
    assert result["reproduction_status"] == "UNEXPECTED"
    assert result["outcome_class"] == "FAILED"
    assert result["invalid_worker_record"]


def test_selected_manifest_is_the_recorded_worker_input(tmp_path):
    selected = tmp_path / "selected-manifest.json"
    selected.write_bytes(reproduction.MANIFEST_PATH.read_bytes())
    output = tmp_path / "worker.json"
    reproduction.worker_main(SimpleNamespace(
        stage="doctor", manifest=selected, output=output,
    ))
    record = json.loads(output.read_text(encoding="utf-8"))
    assert str(selected.resolve()) in record["inputs"]["files"]
    assert str(reproduction.MANIFEST_PATH.resolve()) not in record["inputs"]["files"]


def test_m3a_requires_complete_nonempty_basis_image_evidence():
    expected = reproduction.load_manifest()["stages"]["m3a"]["expected"]
    variants = []
    for name, dimension in expected["closure_dimensions"].items():
        basis = [{"external_pattern": [f"X{index}"], "rop_function": f"q{index}"}
                 for index in range(dimension)]
        images = [{"lhs": row["external_pattern"], "mapped_observable": row["rop_function"],
                   "fresh_fragger_equals_population_dpo_composition": True,
                   "fresh_fragger_equals_molecular_rop_composition_generator": True}
                  for row in basis]
        variants.append({"variant": name, "closure": {"status": "CLOSED",
                         "basis_dimension_including_constant": dimension,
                         "basis_mapping": basis}, "comparable_generator_images": images})
    record = {"variants": variants}
    assert reproduction.validate_m3a(record, expected)[1]["status"] == "PASS"
    record["variants"][0]["comparable_generator_images"] = []
    assert reproduction.validate_m3a(record, expected)[1]["status"] == "FAIL"


def _m3b_record(expected):
    seeds = expected["required_seeds"]
    external, runs = [], []
    for budget in expected["budgets"]:
        rows = [{"status": "EQUAL", "mapped_lhs": {"text": seed},
                 "all_rate_images_equal": True,
                 "external_by_channel": {"k": {"text": seed}},
                 "composition_by_channel": {"k": {"text": seed}}}
                for seed in seeds]
        rows.append({"status": "UNMAPPED"})
        external.append({"status": "PASS", "json": {"seed_names": seeds}})
        runs.append({"budget": budget, "status": "COMPARISON_INCOMPLETE",
                     "mapped_equal": 7, "mapped_mismatch": 0, "unmapped": 1,
                     "total": 8, "rows": rows})
    return {"external_worklists": external, "finite_validation": {"status": "PASS"},
            "symbolic_comparison": {"status": "COMPARISON_INCOMPLETE", "runs": runs}}


def test_m3b_requires_each_seed_equation_and_channel_evidence():
    expected = reproduction.load_manifest()["stages"]["m3b"]["expected"]
    record = _m3b_record(expected)
    assert reproduction.validate_m3b(record, expected)[1]["status"] == "PASS"
    record["symbolic_comparison"]["runs"][0]["rows"] = []
    assert reproduction.validate_m3b(record, expected)[1]["status"] == "FAIL"
    record = _m3b_record(expected)
    record["symbolic_comparison"]["runs"][0]["rows"][0]["mapped_lhs"]["text"] = "replacement"
    assert reproduction.validate_m3b(record, expected)[1]["status"] == "FAIL"


def test_m4_requires_complete_closed_pjc_rows():
    expected = reproduction.load_manifest()["stages"]["m4"]["expected"]
    rows = []
    for key, signature in expected["row_signatures"].items():
        case, gate = key.split(":")
        rows.append({"case": case, "gate": gate == "true", "repetition": 1,
                     "derived_event_count": signature[0], "pjc_differential_order": signature[1],
                     "linear_status": signature[2], "linear_reached_dimension": signature[3],
                     "pjc_status": signature[4], "run_status": signature[5],
                     "oracle_status": signature[6]})
    record = {"status": "COMPLETE", "rows": rows}
    assert reproduction.validate_m4(record, expected, 1)[1]["status"] == "PASS"
    for row in rows:
        row["pjc_status"] = "UNKNOWN"
    assert reproduction.validate_m4(record, expected, 1)[1]["status"] == "FAIL"


def test_tables_use_only_supplied_fresh_records(tmp_path):
    record = {
        "stage": "rop-baseline", "execution_status": "COMPLETE",
        "reproduction_status": "EXPECTED", "elapsed_seconds": 1.0,
        "record_path": str(tmp_path / "stages/rop-baseline/record.json"),
        "analysis_row": {"stage": "rop-baseline", "linear_status": "CLOSED",
                         "linear_dimension": 3, "linear_scope": "fresh",
                         "pjc_status": "CLOSED", "event_count": 2,
                         "operator_order": 1, "pjc_scope": "fresh",
                         "finite_status": "PASS"},
    }
    reproduction.write_tables(tmp_path, [record, {
        "stage": "legacy", "execution_status": "COMPLETE",
        "reproduction_status": "EXPECTED", "record_path": "legacy.json",
        "analysis_row": record["analysis_row"],
    }])
    rows = json.loads((tmp_path / "analysis-rows.json").read_text(encoding="utf-8"))
    assert [row["stage"] for row in rows] == ["rop-baseline"]


def test_tables_keep_failed_analysis_rows_visible(tmp_path):
    reproduction.write_tables(tmp_path, [{
        "stage": "ring8", "execution_status": "PROCESS_ERROR",
        "reproduction_status": "UNEXPECTED", "elapsed_seconds": 0.1,
        "record_path": str(tmp_path / "stages/ring8/record.json"),
        "error": "simulated",
    }])
    rows = json.loads((tmp_path / "analysis-rows.json").read_text(encoding="utf-8"))
    assert rows[0]["stage"] == "ring8"
    assert rows[0]["execution_status"] == "PROCESS_ERROR"
    assert rows[0]["linear_status"] is None


def test_profile_refuses_nonempty_output_directory(tmp_path):
    output = tmp_path / "used"
    output.mkdir()
    (output / "old.json").write_text("{}", encoding="utf-8")
    code = reproduction.run_profile(SimpleNamespace(
        manifest=reproduction.MANIFEST_PATH, output=output, profile="smoke",
        m4_repetitions=1, fragger_checkout=None, java=None,
    ))
    assert code == 2
    assert (output / "old.json").is_file()


def test_stage_exception_is_retained_and_later_stage_runs(tmp_path):
    manifest = reproduction.load_manifest()
    manifest["profiles"]["smoke"] = ["bad", "doctor"]
    manifest["stages"]["bad"] = {"kind": "unsupported", "expected": {"status": "PASS"}}
    selected = tmp_path / "manifest.json"
    selected.write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "out"
    code = reproduction.run_profile(SimpleNamespace(
        manifest=selected, output=output, profile="smoke", m4_repetitions=1,
        m4_timeout=None, fragger_checkout=None, java=None,
    ))
    assert code == 2
    summary = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert [row["stage"] for row in summary["stages"]] == ["bad", "doctor"]
    assert summary["stages"][0]["execution_status"] == "ORCHESTRATION_ERROR"
    assert summary["stages"][1]["reproduction_status"] == "EXPECTED"
