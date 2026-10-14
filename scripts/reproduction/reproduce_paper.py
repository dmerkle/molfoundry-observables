#!/usr/bin/env python3
"""Reproduce the selected analysis evidence with fresh, bounded processes."""
from __future__ import annotations

import mfobs.checks.molecular_rop as molecular_checks
import mfobs.checks.rop_projection as projection_checks

from mfobs.checks import remy_trees as remy_check

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import traceback
from typing import Mapping


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
MANIFEST_PATH = HERE / "paper_evidence_manifest.json"
sys.path.insert(0, str(PROJECT))
sys.path.insert(0, str(PROJECT / "src"))
from mfobs.core.source_identity import source_identity as verified_source_identity
SOURCE_FILES = (
    "scripts/reproduction/paper_evidence_manifest.json",
    "scripts/reproduction/reproduce_paper.py",
    "scripts/reproduction/paper_workers.py",
    "scripts/run_eight_slots.py", "scripts/check_eight_slots.py", "scripts/check_local_slots.py",
    "examples/tcs/compare_conditions.py",
    "examples/tcs/inputs.py",
    "src/mfobs/models/molecular_rop/analysis.py",
    "src/mfobs/models/molecular_rop/counting.py",
    "tests/test_reproduction.py",
    "scripts/benchmarks/rop_scaling.py",
    "scripts/comparison/run_rop_population_comparison.py",
    "scripts/comparison/run_rop_sequence_comparison.py",
    "Dockerfile", "requirements.txt",
    "src/mfobs/example_cli.py",
    "src/mfobs/formatting.py",
    "tests/test_result_formatting.py",
    "examples/chlorination/run.py",
    "examples/growing_ring/run.py",
    "examples/remy/run.py",
    "examples/ring8/run.py",
    "examples/aba/run.py",
    "examples/tcs/replay.py",
    "examples/tcs/_audit.py",
    "examples/tcs/aldol.gml",
    "examples/tcs/query_conditioned.gml",
    "examples/tcs/query_unconditioned.gml",
    "examples/aba/query.json",
    "src/mfobs/core/composition_jump_model.py",
    "src/mfobs/core/finite_language_pjc.py",
    "src/mfobs/core/finite_rule_projection.py",
    "src/mfobs/core/polynomial_jump_closure.py",
    "src/mfobs/core/pattern_support.py",
    "src/mfobs/core/source_identity.py",
    "src/mfobs/core/worklist.py",
    "src/mfobs/core/linear_closure.py",
    "src/mfobs/models/dna_walker/counting.py",
    "src/mfobs/checks/dna_walker.py",
    "src/mfobs/models/chlorination/analysis.py",
    "src/mfobs/models/chlorination/counting.py",
    "src/mfobs/models/chlorination/inputs.py",
    "src/mfobs/models/molecular_rop/inputs.py",
    "src/mfobs/models/aba/analysis.py",
    "src/mfobs/models/aba/counting.py",
    "src/mfobs/models/remy/counting.py",
    "src/mfobs/models/ring8_common.py",
    "src/mfobs/models/structural_nac_common.py",
)
SOURCE_FILES = tuple(dict.fromkeys((*SOURCE_FILES, *(
    path.relative_to(PROJECT).as_posix()
    for folder in ("src/mfobs/models/aba", "src/mfobs/models/molecular_rop",
                   "src/mfobs/models/chlorination", "src/mfobs/models/remy",
                   "src/mfobs/models/dna_walker", "src/mfobs/checks")
    for path in sorted((PROJECT / folder).glob("*.py"))
))))

ANALYSIS_STAGES = {
    "aba",
    "rop-baseline", "rop-noedge", "rop-copolymer", "rop-combined",
    "ring8", "remy",
    "growing-ring",
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_digest(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_manifest(path: Path = MANIFEST_PATH) -> dict:
    value = read_json(path)
    if value.get("schema_version") != 1 or value.get("kind") != "mfobs-evidence-profile":
        raise ValueError("unsupported reproduction manifest")
    stages = value.get("stages")
    if not isinstance(stages, dict) or not stages:
        raise ValueError("manifest must declare stages")
    for profile, names in value.get("profiles", {}).items():
        missing = set(names) - set(stages)
        if missing:
            raise ValueError(f"profile {profile} names missing stages: {sorted(missing)}")
    return value


def git_identity() -> dict:
    return verified_source_identity(PROJECT)


def runtime_identity() -> dict:
    packages = {}
    for name in ("molfoundry", "sympy", "networkx", "numpy", "scipy", "rdkit", "pytest"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {"python": sys.version, "executable": sys.executable,
            "platform": platform.platform(), "machine": platform.machine(),
            "packages": packages}


def source_identity() -> dict:
    files = {name: sha256(PROJECT / name) for name in SOURCE_FILES if (PROJECT / name).is_file()}
    return {"files": files, "sha256": json_digest(files)}


def input_identity(paths: list[Path]) -> dict:
    missing = [str(path.resolve()) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"required input evidence is missing: {missing}")
    files = {str(path.resolve()): sha256(path) for path in paths}
    return {"files": files, "sha256": json_digest(files)}


def selected_manifest(spec: Mapping) -> Path:
    path = Path(spec["_manifest_path"])
    if not path.is_file():
        raise FileNotFoundError(f"selected reproduction manifest is missing: {path}")
    return path


def pjc_event_count(pjc: Mapping) -> int | None:
    if pjc.get("event_class_count") is not None:
        return int(pjc["event_class_count"])
    classes = pjc.get("event_classes")
    return len(classes) if isinstance(classes, list) else None


def pjc_order(pjc: Mapping) -> int | None:
    operator = pjc.get("operator")
    if isinstance(operator, Mapping):
        value = operator.get("maximum_differential_order")
        return int(value) if value is not None else None
    return None


def analysis_summary(stage: str, analysis: Mapping, finite: Mapping | None = None) -> dict:
    linear, pjc = analysis["linear"], analysis["pjc"]
    result = {
        "stage": stage,
        "linear_status": linear["status"],
        "linear_dimension": int(linear.get("dimension", len(linear.get("basis", ())))),
        "linear_scope": linear.get("scope"),
        "pjc_status": pjc["status"],
        "event_count": pjc_event_count(pjc),
        "operator_order": pjc_order(pjc),
        "pjc_scope": pjc.get("scope"),
        "contract_sha256": analysis.get("contract", {}).get("definition_sha256"),
        "composition_status": analysis.get("composition", {}).get("status"),
        "finite_status": finite.get("status") if finite else None,
        "assumption_boundary": analysis.get("contract", {}).get("state_language", {}).get(
            "finite_enumeration_role"
        ),
    }
    return result


def worker_doctor(stage: str, spec: Mapping) -> dict:
    from mfobs.cli import environment
    value = environment()
    status = "PASS" if value.get("engine_loaded") else "FAIL"
    return {"analysis_row": None, "result": value, "observed": {"semantic_status": status},
            "inputs": input_identity([selected_manifest(spec)])}


def worker_tcs_aldol(stage: str, spec: Mapping) -> dict:
    from examples.tcs.replay import run

    result = run()
    counts = result["host_witness"]["input_match_counts"]
    observed = {
        "status": result["status"],
        "conditioned_coefficients": [
            result["variants"]["conditioned"]["selected_net_commutator_coefficients"][side]
            for side in ("before", "after")
        ],
        "conditioned_host_counts": [counts["conditioned"][side]
                                    for side in ("before", "after")],
        "unconditioned_host_counts": [counts["unconditioned"][side]
                                      for side in ("before", "after")],
    }
    paths = [selected_manifest(spec)] + [
        PROJECT / "examples/tcs" / name
        for name in ("replay.py", "_audit.py")
    ]
    paths += sorted((PROJECT / "examples/tcs").glob("*.gml"))
    return {"analysis_row": None, "result": result, "observed": observed,
            "inputs": input_identity(paths)}


def worker_rop(stage: str, spec: Mapping) -> dict:
    from mfobs import ABA, analyze
    from mfobs.models.molecular_rop.inputs import configuration
    selection = spec["selection"]
    cfg = configuration(selection)
    words = cfg.alphabet + tuple("".join(word) for word in cfg.motifs)
    backend = ABA(words, alphabet=cfg.alphabet, gated=cfg.gated)
    analysis = analyze(backend.specification(), backend=backend,
                       linear_limit=int(spec["linear_limit"]), diagnostic_timings=True).record()
    graph = molecular_checks.graph_semantics_validation(selection, int(spec["graph_bound"]))
    projection = molecular_checks.composition_projection_validation(
        selection, int(spec["projection_bound"])
    )
    operator = molecular_checks.polynomial_operator_validation(
        selection, int(spec["operator_bound"])
    )
    reference_digest = projection_checks.analysis_contract(selection).digest
    if any(check.get("input_definition_sha256") != reference_digest
           for check in (graph, projection, operator)):
        raise ValueError("finite reference checks use a different ROP input")
    finite = {
        "status": "PASS" if all(x["status"] == "PASS" for x in (graph, projection, operator)) else "FAIL",
        "graph_semantics": graph, "composition_projection": projection,
        "polynomial_operator": operator,
        "reference_input_definition_sha256": reference_digest,
        "executed_bounds": {"graph_maximum_word_length": spec["graph_bound"],
                            "projection_maximum_word_length": spec["projection_bound"],
                            "operator_maximum_word_length": spec["operator_bound"]},
    }
    row = analysis_summary(stage, analysis, finite)
    model_inputs = [selected_manifest(spec), PROJECT / "src/mfobs/models/molecular_rop/inputs.py"]
    return {"analysis_row": row, "analysis": analysis, "finite_validation": finite,
            "observed": row, "inputs": input_identity(model_inputs)}


def worker_aba(stage: str, spec: Mapping) -> dict:
    """Run the public single-query ROP workflow used by the manuscript."""
    from examples.aba.run import QUERY, run as run_aba

    query = read_json(QUERY)
    requested = {
        "motif": spec["motif"], "alphabet": spec["alphabet"],
        "gate": spec["gated"], "linear_limit": spec["linear_limit"],
        "support_pattern_limit": spec["support_pattern_limit"],
        "support_expansion_limit": spec["support_expansion_limit"],
    }
    if query != requested:
        raise ValueError("ROP paper stage differs from artifact ABA query")
    analysis = run_aba(diagnostic_timings=True)
    support_status = analysis["support_discovery"]["status"]
    from mfobs.checks.terminal_paths import validate_result
    oracle = validate_result(analysis)
    finite = {
        "status": "PASS" if oracle["status"] == "PASS" else "FAIL",
        "support_status": support_status,
        "independent_word_oracle": oracle,
        "finite_enumeration_used_by_analysis": False,
    }
    row = analysis_summary(stage, analysis, finite)
    row["support_status"] = support_status
    model_inputs = [
        selected_manifest(spec),
        QUERY,
        PROJECT / "examples/aba/run.py",
        PROJECT / "src/mfobs/core/pattern_support.py",
        PROJECT / "src/mfobs/models/aba/analysis.py",
        PROJECT / "src/mfobs/models/aba/counting.py",
    ]
    return {"analysis_row": row, "analysis": analysis, "finite_validation": finite,
            "observed": row, "inputs": input_identity(model_inputs)}




def worker_ring8(stage: str, spec: Mapping) -> dict:
    from mfobs.models import ring8_common
    from examples.ring8.run import run_with_evidence
    analysis, fresh = run_with_evidence(worklist_limit=int(spec["worklist_limit"]))
    finite = {
        "status": "PASS" if analysis["finite_rule_projection"]["status"] == "PASS"
                  and analysis["cross_branch_validation"]["status"] == "PASS" else "FAIL",
        "finite_rule_projection": analysis["finite_rule_projection"],
        "cross_branch_validation": analysis["cross_branch_validation"],
        "fresh_worklist": fresh,
    }
    row = analysis_summary(stage, analysis, finite)
    row["fresh_worklist_limit"] = spec["worklist_limit"]
    return {"analysis_row": row, "analysis": analysis, "finite_validation": finite,
            "observed": row, "inputs": input_identity([selected_manifest(spec), ring8_common.MODEL_PATH])}


def worker_remy(stage: str, spec: Mapping) -> dict:
    from examples.remy.run import run
    analysis = run(linear_limit=int(spec["linear_limit"]))
    finite_tree = remy_check.validate_finite(int(spec["finite_tree_bound"]))
    direct = remy_check.validate_molfoundry_direct(int(spec["direct_molfoundry_bound"]))
    finite = {"status": "PASS" if finite_tree["status"] == direct["status"] == "PASS" else "FAIL",
              "recursive_tree_oracle": finite_tree, "direct_molfoundry": direct,
              "theorem_boundary": analysis["pjc"].get("certificate")}
    row = analysis_summary(stage, analysis, finite)
    row["operator_order"] = pjc_order(analysis["pjc"])
    return {"analysis_row": row, "analysis": analysis, "finite_validation": finite,
            "observed": row,
            "inputs": input_identity([selected_manifest(spec), PROJECT / "src/mfobs/models/remy/counting.py"])}


def worker_growing_ring(stage: str, spec: Mapping) -> dict:
    from mfobs.models import structural_nac_common
    from examples.growing_ring.run import run_with_evidence
    model_path = structural_nac_common.EXAMPLE_ROOT / "model.json"
    limit = int(spec["worklist_limit"])
    analysis, fresh = run_with_evidence(worklist_limit=limit)
    row = analysis_summary(stage, analysis)
    row["witness_status"] = analysis["pjc"]["selected_coordinate_non_lumpability"]["status"]
    row["pending_count"] = len(analysis["linear"]["pending"])
    return {"analysis_row": row, "analysis": analysis, "fresh_worklist": fresh,
            "observed": row, "inputs": input_identity([selected_manifest(spec), model_path])}


WORKERS = {"doctor": worker_doctor, "rop": worker_rop,
           "tcs-aldol": worker_tcs_aldol,
           "aba": worker_aba,
           "ring8": worker_ring8,
           "remy": worker_remy, "growing-ring": worker_growing_ring}


from scripts.reproduction.paper_workers import WORKERS as PAPER_WORKERS
WORKERS.update(PAPER_WORKERS)


def compare_expected(observed: Mapping, expected: Mapping) -> dict:
    comparisons = {key: {"expected": value, "actual": observed.get(key),
                         "equal": observed.get(key) == value}
                   for key, value in expected.items()}
    return {"status": "PASS" if all(row["equal"] for row in comparisons.values()) else "FAIL",
            "fields": comparisons}


def worker_main(arguments) -> None:
    started = time.perf_counter()
    record = {"schema_version": 1, "kind": "mfobs-stage-worker",
              "stage": arguments.stage, "execution_status": "ERROR"}
    try:
        manifest = load_manifest(arguments.manifest)
        spec = manifest["stages"][arguments.stage]
        if spec["kind"] != "worker":
            raise ValueError("worker subcommand requires a worker stage")
        worker_spec = {**spec, "_manifest_path": str(arguments.manifest.resolve()),
                       "_output_dir": str(arguments.output.resolve().parent)}
        payload = WORKERS[spec["worker"]](arguments.stage, worker_spec)
        record.update(payload)
        record["expectation_check"] = compare_expected(payload["observed"], spec["expected"])
        record["execution_status"] = "COMPLETE"
        record["reproduction_status"] = (
            "EXPECTED" if record["expectation_check"]["status"] == "PASS" else "UNEXPECTED"
        )
    except Exception as exc:
        record.update(error_type=type(exc).__name__, error=str(exc), traceback=traceback.format_exc())
    record["worker_elapsed_seconds"] = time.perf_counter() - started
    write_json(arguments.output, record)


def descendants(pid: int) -> set[int]:
    try:
        output = subprocess.run(["ps", "-axo", "pid=,ppid="], text=True,
                                capture_output=True, check=True).stdout
    except Exception:
        return set()
    children: dict[int, set[int]] = {}
    for line in output.splitlines():
        child, parent = (int(value) for value in line.split())
        children.setdefault(parent, set()).add(child)
    found, frontier = set(), [pid]
    while frontier:
        current = frontier.pop()
        for child in children.get(current, ()):
            if child not in found:
                found.add(child); frontier.append(child)
    return found


def terminate_tree(process: subprocess.Popen) -> dict:
    found = descendants(process.pid)
    signalled = []
    for pid in sorted(found, reverse=True):
        try:
            os.killpg(pid, signal.SIGKILL); signalled.append(pid)
        except (ProcessLookupError, PermissionError):
            try:
                os.kill(pid, signal.SIGKILL); signalled.append(pid)
            except (ProcessLookupError, PermissionError):
                pass
    try:
        os.killpg(process.pid, signal.SIGKILL); signalled.append(process.pid)
    except ProcessLookupError:
        pass
    return {"signal": "SIGKILL", "descendants_discovered": sorted(found),
            "pids_signalled": signalled}


def run_process(command: list[str], *, stage: str, output: Path, timeout: float,
                required_output: Path | None = None) -> dict:
    stage_dir = output / "stages" / stage
    stage_dir.mkdir(parents=True, exist_ok=True)
    stdout_path, stderr_path = stage_dir / "stdout.log", stage_dir / "stderr.log"
    print(f"[{stage}] starting (deadline {timeout:.0f}s)", flush=True)
    started = time.perf_counter()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(PROJECT / "src") + (
        os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else ""
    )
    timed_out, termination = False, None
    with stdout_path.open("w", encoding="utf-8") as stdout_stream, stderr_path.open("w", encoding="utf-8") as stderr_stream:
        process = subprocess.Popen(
            command, cwd=PROJECT, text=True, stdout=stdout_stream,
            stderr=stderr_stream, start_new_session=True, env=environment,
        )
        deadline = started + timeout
        while process.poll() is None:
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                timed_out = True; termination = terminate_tree(process)
                process.wait()
                break
            try:
                process.wait(timeout=min(15.0, remaining))
            except subprocess.TimeoutExpired:
                print(f"[{stage}] still running after {time.perf_counter()-started:.1f}s", flush=True)
    elapsed = time.perf_counter() - started
    status = "TIMEOUT" if timed_out else (
        "PROCESS_ERROR" if process.returncode != 0 else
        "INVALID_OUTPUT" if required_output is not None and not required_output.is_file()
        else "COMPLETE"
    )
    record = {"execution_status": status, "exact_command": command,
              "timeout_seconds": timeout, "elapsed_seconds": elapsed,
              "process_returncode": process.returncode, "stdout_log": str(stdout_path),
              "stderr_log": str(stderr_path), "required_output": str(required_output) if required_output else None}
    if termination: record["termination"] = termination
    print(f"[{stage}] {status} in {elapsed:.3f}s", flush=True)
    return record


def validate_m3a(value: Mapping, expected: Mapping) -> tuple[dict, dict]:
    dimensions = {row["variant"]: row["closure"]["basis_dimension_including_constant"]
                  for row in value["variants"]}
    statuses = {row["variant"]: row["closure"]["status"] for row in value["variants"]}
    coverage_ok = True
    equalities = []
    for row in value["variants"]:
        images = row.get("comparable_generator_images")
        basis = row["closure"].get("basis_mapping")
        if not isinstance(images, list) or not images or not isinstance(basis, list):
            coverage_ok = False
            continue
        expected_pairs = [(tuple(item["external_pattern"]), item["rop_function"]) for item in basis]
        image_pairs = [(tuple(item["lhs"]), item["mapped_observable"]) for item in images]
        if (len(images) != row["closure"]["basis_dimension_including_constant"]
                or len(set(image_pairs)) != len(image_pairs)
                or sorted(image_pairs) != sorted(expected_pairs)):
            coverage_ok = False
        for item in images:
            for field in ("fresh_fragger_equals_population_dpo_composition",
                          "fresh_fragger_equals_molecular_rop_composition_generator"):
                equalities.append(item.get(field) is True)
    observed = {"status": "PASS" if coverage_ok and equalities and all(equalities) else "FAIL",
                "closure_dimensions": dimensions, "closure_statuses": statuses}
    return observed, compare_expected(observed, expected)


def validate_m3b(value: Mapping, expected: Mapping) -> tuple[dict, dict]:
    runs = value["symbolic_comparison"]["runs"]
    required_seeds = expected["required_seeds"]
    evidence_ok = len(runs) == len(value["external_worklists"]) == len(expected["budgets"])
    for external, run in zip(value["external_worklists"], runs):
        if external.get("json", {}).get("seed_names") != required_seeds:
            evidence_ok = False
        rows = run.get("rows")
        if not isinstance(rows, list):
            evidence_ok = False
            continue
        equal_rows = [row for row in rows if row.get("status") == "EQUAL"]
        mismatch_rows = [row for row in rows if row.get("status") == "MISMATCH"]
        unmapped_rows = [row for row in rows if row.get("status") == "UNMAPPED"]
        if (len(equal_rows) != run.get("mapped_equal")
                or len(mismatch_rows) != run.get("mapped_mismatch")
                or len(unmapped_rows) != run.get("unmapped")
                or len(rows) != run.get("total")):
            evidence_ok = False
        for seed in required_seeds:
            matches = [row for row in rows
                       if isinstance(row.get("mapped_lhs"), Mapping)
                       and row["mapped_lhs"].get("text") == seed]
            if len(matches) != 1:
                evidence_ok = False
                continue
            row = matches[0]
            external_channels = row.get("external_by_channel")
            composition_channels = row.get("composition_by_channel")
            if (row.get("status") != "EQUAL" or row.get("all_rate_images_equal") is not True
                    or not isinstance(external_channels, dict) or not external_channels
                    or external_channels != composition_channels):
                evidence_ok = False
    observed = {"status": value["symbolic_comparison"]["status"],
                "minimum_seed_equalities": min(row["mapped_equal"] for row in runs),
                "mapped_mismatches": sum(row["mapped_mismatch"] for row in runs),
                "finite_status": value["finite_validation"]["status"],
                "budgets": [row["budget"] for row in runs],
                "external_statuses": [row["status"] for row in value["external_worklists"]],
                "evidence_status": "PASS" if evidence_ok else "FAIL",
                "required_seeds": required_seeds}
    return observed, compare_expected(observed, expected)


def validate_m4(value: Mapping, expected: Mapping, repetitions: int) -> tuple[dict, dict]:
    grouped = {}
    for row in value.get("rows", ()):
        key = f"{row['case']}:{str(row['gate']).lower()}"
        signature = [row.get("derived_event_count"), row.get("pjc_differential_order"),
                     row.get("linear_status"), row.get("linear_reached_dimension"),
                     row.get("pjc_status"), row.get("run_status"), row.get("oracle_status")]
        grouped.setdefault(key, []).append((row.get("repetition"), signature))
    signatures = {
        key: rows[0][1] for key, rows in grouped.items()
        if len(rows) == repetitions
        and {row[0] for row in rows} == set(range(1, repetitions + 1))
        and all(row[1] == rows[0][1] for row in rows)
    }
    valid_rows = signatures == expected["row_signatures"]
    observed = {"status": value.get("status") if valid_rows else "INVALID_EVIDENCE",
                "cases_per_repetition": len(value.get("rows", ())) // repetitions,
                "row_signatures": signatures}
    return observed, compare_expected(observed, expected)


def validate_worker_record(value: Mapping, stage: str) -> None:
    required = {"schema_version", "kind", "stage", "execution_status",
                "worker_elapsed_seconds"}
    missing = required - set(value)
    if missing:
        raise ValueError(f"worker record lacks mandatory fields: {sorted(missing)}")
    if (value["schema_version"] != 1 or value["kind"] != "mfobs-stage-worker"
            or value["stage"] != stage):
        raise ValueError("worker record schema, kind, or stage does not match invocation")
    if value["execution_status"] == "COMPLETE":
        complete_fields = {"reproduction_status", "expectation_check", "observed", "inputs"}
        missing = complete_fields - set(value)
        if missing:
            raise ValueError(f"complete worker record lacks fields: {sorted(missing)}")
    elif value["execution_status"] == "ERROR":
        if not value.get("error") or not value.get("error_type"):
            raise ValueError("failed worker record lacks error evidence")
    else:
        raise ValueError(f"unknown worker execution status: {value['execution_status']!r}")


def _run_stage(stage: str, spec: Mapping, arguments, output: Path) -> dict:
    stage_dir = output / "stages" / stage
    stage_dir.mkdir(parents=True, exist_ok=True)
    record_path = stage_dir / "record.json"
    if spec["kind"] == "worker":
        record_path.unlink(missing_ok=True)
        command = [sys.executable, str(Path(__file__).resolve()), "worker", "--stage", stage,
                   "--manifest", str(arguments.manifest.resolve()), "--output", str(record_path)]
        process = run_process(command, stage=stage, output=output,
                              timeout=float(spec["timeout_seconds"]), required_output=record_path)
        worker = {}
        invalid_worker = None
        if record_path.is_file():
            try:
                worker = read_json(record_path)
                validate_worker_record(worker, stage)
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                invalid_worker = f"{type(exc).__name__}: {exc}"
        result = {**worker, **process, "record_path": str(record_path)}
        if invalid_worker:
            result.update(execution_status="INVALID_OUTPUT", reproduction_status="UNEXPECTED",
                          invalid_worker_record=invalid_worker)
        elif process["execution_status"] != "COMPLETE":
            result["reproduction_status"] = "UNEXPECTED"
        elif worker.get("execution_status") != "COMPLETE":
            result["execution_status"] = worker.get("execution_status", "INVALID_OUTPUT")
            result["reproduction_status"] = "UNEXPECTED"
    elif spec["kind"] == "pytest":
        command = [sys.executable, "-m", "pytest", "-q", "-o",
                   "cache_dir=/tmp/pytest-cache", *spec["tests"]]
        process = run_process(command, stage=stage, output=output,
                              timeout=float(spec["timeout_seconds"]))
        observed = {"status": "PASS" if process["execution_status"] == "COMPLETE" else "FAIL"}
        result = {**process, "observed": observed,
                  "expectation_check": compare_expected(observed, spec["expected"]),
                  "record_path": str(record_path)}
        result["reproduction_status"] = "EXPECTED" if result["expectation_check"]["status"] == "PASS" else "UNEXPECTED"
        result["inputs"] = input_identity([PROJECT / name.split("::", 1)[0] for name in spec["tests"]])
    elif spec["kind"] in {"m3a", "m3b"}:
        missing = [name for name in spec["required_arguments"]
                   if getattr(arguments, name) is None or not getattr(arguments, name).exists()]
        if missing:
            result = {"execution_status": "INCOMPLETE_PREREQUISITE",
                      "reproduction_status": "INCOMPLETE", "missing_arguments": missing,
                      "record_path": str(record_path)}
        else:
            script = "run_rop_population_comparison.py" if spec["kind"] == "m3a" else "run_rop_sequence_comparison.py"
            raw_dir = stage_dir / "raw"
            command = [sys.executable, str(PROJECT / "scripts/comparison" / script),
                       "--checkout", str(arguments.fragger_checkout.resolve()),
                       "--java", str(arguments.java.resolve()), "--out", str(raw_dir)]
            raw = raw_dir / spec["required_output"]
            process = run_process(command, stage=stage, output=output,
                                  timeout=float(spec["timeout_seconds"]), required_output=raw)
            result = {**process, "record_path": str(record_path), "raw_result": str(raw)}
            result["inputs"] = input_identity([
                arguments.fragger_checkout / "compiler" / name
                for name in ("scala-compiler-2.13.2.jar", "scala-library-2.13.2.jar", "scala-reflect-2.13.2.jar")
            ] + [arguments.java, PROJECT / "scripts/comparison" / script])
            if raw.is_file() and process["execution_status"] == "COMPLETE":
                observed, check = (validate_m3a(read_json(raw), spec["expected"])
                                   if spec["kind"] == "m3a" else
                                   validate_m3b(read_json(raw), spec["expected"]))
                result.update(observed=observed, expectation_check=check,
                              reproduction_status="EXPECTED" if check["status"] == "PASS" else "UNEXPECTED")
    elif spec["kind"] == "m4":
        raw_dir = stage_dir / "raw"
        inner_timeout = arguments.m4_timeout or spec["inner_timeout_seconds"]
        command = [sys.executable, str(PROJECT / "scripts/benchmarks/rop_scaling.py"),
                   "run", "--profile", "full", "--repetitions", str(arguments.m4_repetitions),
                   "--timeout", str(inner_timeout), "--output", str(raw_dir)]
        timeout = arguments.m4_repetitions * spec["expected"]["cases_per_repetition"] * inner_timeout + 300
        raw = raw_dir / spec["required_output"]
        process = run_process(command, stage=stage, output=output, timeout=timeout, required_output=raw)
        result = {**process, "record_path": str(record_path), "raw_result": str(raw)}
        result["inputs"] = input_identity([
            arguments.manifest, PROJECT / "scripts/benchmarks/rop_scaling.py",
            PROJECT / "src/mfobs/models/molecular_rop/inputs.py",
        ])
        if raw.is_file() and process["execution_status"] == "COMPLETE":
            observed, check = validate_m4(read_json(raw), spec["expected"], arguments.m4_repetitions)
            result.update(observed=observed, expectation_check=check,
                          reproduction_status="EXPECTED" if check["status"] == "PASS" else "UNEXPECTED")
    else:
        raise ValueError(f"unknown stage kind: {spec['kind']}")
    result.setdefault("reproduction_status", "UNEXPECTED")
    expected_is_incomplete = any(
        isinstance(value, str) and (value.startswith("UNKNOWN_") or value == "COMPARISON_INCOMPLETE")
        for value in spec.get("expected", {}).values()
    )
    result["outcome_class"] = (
        "EXPECTED_INCOMPLETE" if result["reproduction_status"] == "EXPECTED" and expected_is_incomplete
        else "PASS" if result["reproduction_status"] == "EXPECTED" else "FAILED"
    )
    result.update(schema_version=1, kind="mfobs-stage-record", stage=stage,
                  declaration=dict(spec), source_identity=source_identity(),
                  runtime_identity=runtime_identity(), project_revision=git_identity())
    write_json(record_path, result)
    return result


def run_stage(stage: str, spec: Mapping, arguments, output: Path) -> dict:
    """Retain one failed record instead of aborting the declared profile."""
    try:
        return _run_stage(stage, spec, arguments, output)
    except Exception as exc:
        record_path = output / "stages" / stage / "record.json"
        result = {
            "schema_version": 1, "kind": "mfobs-stage-record",
            "stage": stage, "declaration": dict(spec),
            "execution_status": "ORCHESTRATION_ERROR",
            "reproduction_status": "UNEXPECTED", "outcome_class": "FAILED",
            "error_type": type(exc).__name__, "error": str(exc),
            "traceback": traceback.format_exc(), "record_path": str(record_path),
        }
        try:
            result.update(source_identity=source_identity(), runtime_identity=runtime_identity(),
                          project_revision=git_identity())
        except Exception as identity_error:
            result["identity_error"] = f"{type(identity_error).__name__}: {identity_error}"
        write_json(record_path, result)
        return result


def relative_link(path: str, report: Path) -> str:
    try:
        return os.path.relpath(path, report.parent)
    except ValueError:
        return path


def write_tables(output: Path, records: list[dict]) -> None:
    rows = []
    for record in records:
        if record["stage"] not in ANALYSIS_STAGES and (
                not record.get("analysis_row") or record["analysis_row"].get("stage") != record["stage"]):
            continue
        row = {"stage": record["stage"], "linear_status": None,
               "linear_dimension": None, "linear_scope": None,
               "pjc_status": None, "event_count": None, "operator_order": None,
               "pjc_scope": None, "finite_status": None}
        row.update(record.get("analysis_row") or {})
        row["execution_status"] = record["execution_status"]
        row["reproduction_status"] = record["reproduction_status"]
        row["error"] = record.get("error") or record.get("invalid_worker_record")
        row["elapsed_seconds"] = record.get("elapsed_seconds")
        row["raw_result"] = record["record_path"]
        rows.append(row)
    write_json(output / "analysis-rows.json", rows)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (output / "analysis-rows.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    report = output / "evidence-table.md"
    lines = ["# Fresh paper-evidence analysis rows", "",
             "Expected `UNKNOWN_*` values reproduce bounded outcomes; they do not mean closure.", "",
             "| Model | Linear status (dimension) | Linear scope | pJC status | Events/order | pJC scope | Finite check | Raw |",
             "|---|---|---|---|---|---|---|---|"]
    for row in rows:
        raw = relative_link(row["raw_result"], report)
        lines.append("| " + " | ".join([
            row["stage"], (f"{row['linear_status']} ({row['linear_dimension']})"
                           if row.get("linear_status") else row["execution_status"]),
            str(row.get("linear_scope") or "—"), str(row.get("pjc_status") or "—"),
            f"{row.get('event_count') if row.get('event_count') is not None else '—'}/{row.get('operator_order') if row.get('operator_order') is not None else '—'}",
            str(row.get("pjc_scope") or "—"), str(row.get("finite_status") or "—"),
            f"[JSON]({raw})",
        ]) + " |")
    lines += ["", "Rémy's finite checks are validation only; its all-tree P2/P3 result uses the supplied theorem boundary. Growing-ring UNKNOWN retains its pending frontier and a separate non-lumpability witness.", ""]
    report.write_text("\n".join(lines), encoding="utf-8")
    claims = ["# Checks against the tagged paper", "",
              "Each row links to fresh output, its expected values and all observed values.", "",
              "| Stage | Outcome | Observed checks |", "|---|---|---|"]
    for record in records:
        raw = relative_link(record["record_path"], output / "claims.md")
        observed = json.dumps(record.get("observed", {}), sort_keys=True).replace("|", "&#124;")
        claims.append(f"| [{record['stage']}]({raw}) | {record.get('outcome_class', record['reproduction_status'])} | `{observed}` |")
    (output / "claims.md").write_text("\n".join(claims) + "\n", encoding="utf-8")


def write_inventory(output: Path) -> None:
    files = {}
    inventory = output / "sha256-inventory.json"
    for path in sorted(output.rglob("*")):
        if path.is_file() and path != inventory:
            files[path.relative_to(output).as_posix()] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    write_json(inventory, {"schema_version": 1, "scope": "all retained output files except this self-referential inventory", "files": files})


def run_profile(arguments) -> int:
    started = time.perf_counter()
    manifest = load_manifest(arguments.manifest)
    output = arguments.output.resolve()
    if output.exists() and any(output.iterdir()):
        print(f"output directory must be new or empty: {output}", file=sys.stderr)
        return 2
    output.mkdir(parents=True, exist_ok=True)
    initial = {"schema_version": 1, "kind": "mfobs-reproduction-run",
               "profile": arguments.profile, "manifest": str(arguments.manifest.resolve()),
               "manifest_sha256": sha256(arguments.manifest),
               "stage_order": manifest["profiles"][arguments.profile],
               "coverage_gaps": manifest.get("coverage_gaps", {}).get(arguments.profile, []),
               "arguments": {"m4_repetitions": arguments.m4_repetitions,
                             "m4_timeout": arguments.m4_timeout,
                             "fragger_checkout": str(arguments.fragger_checkout) if arguments.fragger_checkout else None,
                             "java": str(arguments.java) if arguments.java else None},
               "project_revision_before": git_identity(), "source_identity": source_identity(),
               "runtime_identity": runtime_identity()}
    write_json(output / "run.json", initial)
    records = []
    for stage in initial["stage_order"]:
        records.append(run_stage(stage, manifest["stages"][stage], arguments, output))
    counts = {"passed": sum(row["outcome_class"] == "PASS" for row in records),
              "expected_incomplete": sum(row["outcome_class"] == "EXPECTED_INCOMPLETE" for row in records),
              "failed": sum(row["outcome_class"] == "FAILED" for row in records)}
    complete = counts["failed"] == 0
    initial.update(
        stages=[{"stage": row["stage"], "execution_status": row["execution_status"],
                 "reproduction_status": row["reproduction_status"],
                 "outcome_class": row["outcome_class"],
                 "record_path": row["record_path"], "elapsed_seconds": row.get("elapsed_seconds")}
                for row in records],
        counts=counts, elapsed_seconds=time.perf_counter() - started,
        project_revision_after=git_identity(),
        status="REPRODUCED_EXPECTED_OUTCOMES" if complete else "INCOMPLETE_OR_UNEXPECTED",
        interpretation=("Expected UNKNOWN and declared partial coverage are included only when their structured assertions pass; this status does not assert universal closure."),
    )
    queries = [r for r in records if r["stage"].startswith("query-") and not r["declaration"].get("gated")]
    rule_hashes = [r.get("observed", {}).get("source_rules_sha256") for r in queries]
    initial["query_reuse_same_source_rules"] = (
        bool(rule_hashes) and None not in rule_hashes and len(set(rule_hashes)) == 1
    ) if queries else None
    if queries and not initial["query_reuse_same_source_rules"]:
        complete = False
        initial["status"] = "INCOMPLETE_OR_UNEXPECTED"
    initial["paper_claim_coverage"] = "PARTIAL" if initial["coverage_gaps"] else "DECLARED_PROFILE"
    write_tables(output, records)
    write_json(output / "source-runtime-input-manifest.json", {
        "schema_version": 1, "project_revision": initial["project_revision_after"],
        "source": source_identity(), "runtime": runtime_identity(),
        "manifest": {"path": str(arguments.manifest.resolve()), "sha256": sha256(arguments.manifest)},
        "stage_inputs": {row["stage"]: row.get("inputs") for row in records},
    })
    write_json(output / "run.json", initial); write_inventory(output)
    print(f"profile: {arguments.profile} -> {initial['status']} in {initial['elapsed_seconds']:.3f}s", flush=True)
    print(f"summary: {output / 'run.json'}", flush=True)
    for gap in initial["coverage_gaps"]:
        print(f"COVERAGE GAP: {gap['claim']}: {gap['reason']}", flush=True)
    return 0 if complete else 2


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    worker = commands.add_parser("worker", help=argparse.SUPPRESS)
    worker.add_argument("--stage", required=True); worker.add_argument("--manifest", type=Path, required=True)
    worker.add_argument("--output", type=Path, required=True)
    run = commands.add_parser("run")
    run.add_argument("--profile", choices=("smoke", "paper", "supplementary"), required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    run.add_argument("--fragger-checkout", type=Path); run.add_argument("--java", type=Path)
    run.add_argument("--m4-repetitions", type=int, default=1)
    run.add_argument("--m4-timeout", type=float,
                     help="per-process M4 wall deadline; scientific limits are unchanged")
    return result


def main() -> None:
    arguments = parser().parse_args()
    if arguments.command == "worker": worker_main(arguments); return
    if arguments.m4_repetitions < 1: raise SystemExit("--m4-repetitions must be positive")
    raise SystemExit(run_profile(arguments))


if __name__ == "__main__":
    main()
