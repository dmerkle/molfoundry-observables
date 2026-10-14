"""Small protocol tests for the isolated ROP scaling runner."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace


sys.path.insert(0, str(Path(__file__).parent))
from scripts.benchmarks import rop_scaling


def test_full_matrix_is_nested_and_has_eight_fixed_cases():
    rows = rop_scaling.schedule(
        profile="full", case=None, gate="both", repetitions=1
    )
    assert len(rows) == 8
    assert [(row["case"], row["gate"]) for row in rows] == [
        (case, gate) for case in rop_scaling.CASES for gate in (False, True)
    ]
    motifs = [set(rop_scaling.CASES[name]["motifs"]) for name in rop_scaling.CASES]
    assert all(left <= right for left, right in zip(motifs, motifs[1:]))
    assert rop_scaling.CASES["abba"] == {
        "motifs": ("AA", "AB", "BA", "BB", "ABA", "ABBA"),
        "context_depth": 3,
    }
    assert [
        rop_scaling.EXPECTED_EVENT_COUNTS[(row["case"], row["gate"])]
        for row in rows
    ] == [4, 6, 8, 10, 16, 18, 32, 34]


def test_smoke_and_case_selection_remain_explicit():
    smoke = rop_scaling.schedule(
        profile="smoke", case=None, gate="both", repetitions=1
    )
    assert [(row["case"], row["gate"]) for row in smoke] == [
        ("letters", False),
        ("letters", True),
    ]
    selected = rop_scaling.schedule(
        profile="full", case="aba", gate="true", repetitions=2
    )
    assert [(row["case"], row["gate"], row["repetition"]) for row in selected] == [
        ("aba", True, 1),
        ("aba", True, 2),
    ]


def test_letters_worker_keeps_complete_suffix_workload(tmp_path):
    output = tmp_path / "letters.json"
    rop_scaling.worker(SimpleNamespace(
        case="letters", gate=False, repetition=1, output=output,
    ))
    record = json.loads(output.read_text(encoding="utf-8"))
    assert record["run_status"] == "COMPLETE", record.get("traceback")
    events = record["analysis"]["composition"]["event_classes"]
    assert len(events) == 4
    assert all("activity_polynomial" in event for event in events)


def complete_record() -> dict:
    return {
        "case": "letters",
        "motifs": [],
        "gate": False,
        "context_depth": 0,
        "repetition": 1,
        "run_status": "COMPLETE",
        "linear_limit": 24,
        "record_path": "/tmp/one.json",
        "total_process_wall_seconds": 1.0,
        "analysis": {
            "configuration": {"sha256": "config"},
            "contract": {"definition_sha256": "contract"},
            "composition": {
                "event_classes": [
                    {"activity_polynomial": {"degree": 0}},
                    {"activity_polynomial": {"degree": 2}},
                ]
            },
            "linear": {
                "status": "UNKNOWN_DIMENSION_LIMIT",
                "basis": [f"O{index}" for index in range(24)],
                "origins": {
                    "N_A": {"origin": "selected-seed"},
                    "M_A": {"origin": "generator-image-residual"},
                },
            },
            "pjc": {
                "status": "CLOSED",
                "observable_coordinates": ["N_A", "M_A"],
                "operator": {"maximum_differential_order": 2},
            },
            "word_oracle_comparison": {"status": "PASS"},
            "diagnostics": {
                "seconds": {
                    "total": 0.9,
                    "configuration_and_contract": 0.1,
                    "composition": 0.4,
                    "independent_oracle_comparison": 0.1,
                    "linear_discovery": 0.2,
                    "polynomial_jump_closure": 0.1,
                }
            },
        },
    }


def test_result_extraction_uses_analysis_fields_not_case_expectations():
    row = rop_scaling.extract_row(complete_record())
    assert row["selected_seed_count"] == 1
    assert row["coordinate_count"] == 2
    assert row["derived_event_count"] == 2
    assert row["maximum_propensity_degree"] == 2
    assert row["pjc_differential_order"] == 2
    assert row["linear_status"] == "UNKNOWN_DIMENSION_LIMIT"
    assert row["linear_reached_dimension"] == 24
    assert row["oracle_status"] == "PASS"


def test_summary_retains_timeout_and_error_rows():
    complete = complete_record()
    failure_spec = {
        "case": "dyads",
        "motifs": ["AA", "AB", "BA", "BB"],
        "gate": False,
        "context_depth": 1,
        "repetition": 1,
    }
    timeout = rop_scaling.failure_record(
        failure_spec, "TIMEOUT", "simulated"
    )
    error = rop_scaling.failure_record(
        {**failure_spec, "case": "aba", "context_depth": 2},
        "ERROR",
        "simulated",
    )
    protocol = {"requested_run_count": 3, "requested_repetitions": 1}
    summary = rop_scaling.summary_record([complete, timeout, error], protocol)
    assert summary["status"] == "INCOMPLETE"
    assert summary["recorded_run_count"] == 3
    assert [row["run_status"] for row in summary["rows"]] == [
        "COMPLETE",
        "TIMEOUT",
        "ERROR",
    ]
