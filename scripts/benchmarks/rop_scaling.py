#!/usr/bin/env python3
"""Fresh-process scaling measurements for nested ROP observable queries.

The parent runs workers serially and enforces a real process deadline.  Each
worker constructs one configuration, compiles its complete suffix refinement,
and runs the shared linear and jump analysis with linear budget 24. This fixed
benchmark workload is separate from query-driven support discovery.
"""
from __future__ import annotations

import mfobs.checks.rop_projection as rop_rop_projection

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

from mfobs.core.source_identity import source_identity as verified_source_identity


PROCESS_STARTED = time.perf_counter()
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
LINEAR_LIMIT = 24
DEFAULT_TIMEOUT_SECONDS = 600.0
CASES = {
    "letters": {"motifs": (), "context_depth": 0},
    "dyads": {"motifs": ("AA", "AB", "BA", "BB"), "context_depth": 1},
    "aba": {
        "motifs": ("AA", "AB", "BA", "BB", "ABA"),
        "context_depth": 2,
    },
    "abba": {
        "motifs": ("AA", "AB", "BA", "BB", "ABA", "ABBA"),
        "context_depth": 3,
    },
}
EXPECTED_EVENT_COUNTS = {
    ("letters", False): 4,
    ("letters", True): 6,
    ("dyads", False): 8,
    ("dyads", True): 10,
    ("aba", False): 16,
    ("aba", True): 18,
    ("abba", False): 32,
    ("abba", True): 34,
}
SOURCE_FILES = (
    "src/mfobs/models/molecular_rop/inputs.py",
    "src/mfobs/core/molecules.py",
    "src/mfobs/core/composition_jump_model.py",
    "src/mfobs/core/polynomial_jump_closure.py",
    "src/mfobs/core/counting_linear_closure.py",
    "src/mfobs/core/linear_closure.py",
    "src/mfobs/models/aba/inputs.py",
    "src/mfobs/models/aba/counting.py",
    "src/mfobs/models/aba/projection.py",
    "src/mfobs/core/engine.py",
    "scripts/benchmarks/rop_scaling.py",
)


def digest_json(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def source_hashes() -> dict:
    files = {
        name: hashlib.sha256((PROJECT / name).read_bytes()).hexdigest()
        for name in SOURCE_FILES
    }
    return {"sha256": digest_json(files), "files": files}


def git_revision() -> dict:
    return verified_source_identity(PROJECT)


def runtime_identity() -> dict:
    distribution = importlib.metadata.distribution("molfoundry")
    files = {}
    for relative in distribution.files or ():
        normalized = str(relative).replace("\\", "/")
        path = Path(distribution.locate_file(relative))
        if path.is_file() and (
            path.suffix.lower() in {".so", ".dylib", ".pyd", ".dll", ".js"}
            or normalized == "molfoundry/_engine.py"
        ):
            files[normalized] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "distribution": "molfoundry",
        "version": distribution.version,
        "build_sha256": digest_json(files),
        "build_files": files,
    }


def host_environment() -> dict:
    uname = platform.uname()
    return {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "system": uname.system,
        "release": uname.release,
        "machine": uname.machine,
        "processor": uname.processor,
        "logical_cpu_count": os.cpu_count(),
    }


def selected_cases(profile: str, case: str | None) -> tuple[str, ...]:
    if case is not None:
        return (case,)
    if profile == "smoke":
        return ("letters",)
    return tuple(CASES)


def schedule(
    *, profile: str, case: str | None, gate: str, repetitions: int
) -> list[dict]:
    if repetitions < 1:
        raise ValueError("repetitions must be positive")
    gates = {"false": (False,), "true": (True,), "both": (False, True)}[gate]
    return [
        {
            "case": case_name,
            "motifs": list(CASES[case_name]["motifs"]),
            "context_depth": CASES[case_name]["context_depth"],
            "gate": gated,
            "repetition": repetition,
        }
        for repetition in range(1, repetitions + 1)
        for case_name in selected_cases(profile, case)
        for gated in gates
    ]


def measure_refined_query(selection):
    """Keep the original complete-suffix workload for timing comparisons."""
    from mfobs.models.molecular_rop import inputs as rop
    from mfobs.core.composition_jump_model import compare_polynomial_jump_models
    from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
    from mfobs.core.counting_linear_closure import discover_linear_observables
    from mfobs.models.aba import inputs as paths
    from mfobs.models.aba.counting import PathCounter
    from mfobs.models.aba.projection import TerminalWordGraphModel

    started, phase = time.perf_counter(), time.perf_counter()
    effective = rop.effective_configuration(selection, linear_limit=LINEAR_LIMIT)
    contract, language = rop_rop_projection.analysis_contract(selection), rop_rop_projection.grammar(selection)
    timings = {"configuration_and_contract": time.perf_counter() - phase}
    if TerminalWordGraphModel(language).refinement_consistency_check()["status"] != "PASS":
        raise AssertionError("inconsistent terminal-word refinement")
    phase = time.perf_counter()
    queries = {p.name: p for p in map(paths.query_pattern,
        [*selection.alphabet, *("".join(word) for word in selection.motifs)])}
    linear = discover_linear_observables(
        paths.source_rules(selection.alphabet, selection.gated), queries,
        counter=PathCounter(state_space=paths.state_family(selection.alphabet, selection.gated)),
        maximum_dimension=LINEAR_LIMIT, input_definition_sha256=contract.digest)
    timings["linear_discovery"] = time.perf_counter() - phase
    phase = time.perf_counter()
    compilation = rop_rop_projection.composition_compilation(selection)
    timings["composition"] = time.perf_counter() - phase
    phase = time.perf_counter()
    oracle = compare_polynomial_jump_models(compilation.model, language.polynomial_jump_model())
    timings["independent_oracle_comparison"] = time.perf_counter() - phase
    if oracle["status"] != "PASS":
        raise AssertionError("composition-derived model differs from the word oracle")
    phase = time.perf_counter()
    pjc = analyze_polynomial_jump_closure(compilation.model)
    timings["polynomial_jump_closure"] = time.perf_counter() - phase
    timings["total"] = time.perf_counter() - started
    return {"configuration": effective, "contract": contract.record(),
            "composition": compilation.evidence, "word_oracle_comparison": oracle,
            "linear": linear, "pjc": pjc, "diagnostics": {"seconds": timings}}


def worker(arguments) -> None:
    record = {
        "schema_version": 1,
        "kind": "rop-scaling-worker-record",
        "case": arguments.case,
        "motifs": list(CASES[arguments.case]["motifs"]),
        "context_depth": CASES[arguments.case]["context_depth"],
        "gate": arguments.gate,
        "repetition": arguments.repetition,
        "linear_limit": LINEAR_LIMIT,
        "pid": os.getpid(),
        "analysis_call_count": 0,
        "host": host_environment(),
    }
    try:
        import_started = time.perf_counter()
        from mfobs.models.molecular_rop import inputs as copolymer_rop

        record["worker_startup_and_import_seconds"] = (
            time.perf_counter() - PROCESS_STARTED
        )
        record["model_import_seconds"] = time.perf_counter() - import_started
        record["runtime_identity"] = runtime_identity()

        selection_started = time.perf_counter()
        selection = copolymer_rop.configured_variant(
            motifs=CASES[arguments.case]["motifs"],
            gated=arguments.gate,
            alphabet=("A", "B"),
            title=f"ROP scaling · {arguments.case}",
        )
        record["selection_construction_seconds"] = (
            time.perf_counter() - selection_started
        )

        analysis_started = time.perf_counter()
        record["analysis_call_count"] = 1
        analysis = measure_refined_query(selection)
        record["analysis_call_wall_seconds"] = time.perf_counter() - analysis_started
        record["analysis"] = analysis
        derived_events = len(analysis["composition"]["event_classes"])
        expected_events = EXPECTED_EVENT_COUNTS[(arguments.case, arguments.gate)]
        record["event_count_review_sanity"] = {
            "status": "PASS" if derived_events == expected_events else "FAIL",
            "expected": expected_events,
            "actual": derived_events,
            "role": "post-analysis review sanity check; never supplied to computation",
        }
        record["run_status"] = "COMPLETE"
    except Exception as exc:
        record.update(
            run_status="ERROR",
            error_type=type(exc).__name__,
            error=str(exc),
            traceback=traceback.format_exc(),
        )
    finally:
        record["worker_wall_seconds"] = time.perf_counter() - PROCESS_STARTED
        write_json(arguments.output, record)


def failure_record(spec: dict, status: str, message: str) -> dict:
    return {
        "schema_version": 1,
        "kind": "rop-scaling-worker-record",
        **spec,
        "linear_limit": LINEAR_LIMIT,
        "run_status": status,
        "error": message,
    }


def run_process(
    spec: dict, *, index: int, total: int, output: Path, timeout: float
) -> dict:
    stem = (
        f"{index:02d}-{spec['case']}-gate-{str(spec['gate']).lower()}"
        f"-r{spec['repetition']}"
    )
    record_path = output / f"{stem}.json"
    stdout_path = output / f"{stem}.stdout.log"
    stderr_path = output / f"{stem}.stderr.log"
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "worker",
        "--case",
        spec["case"],
        "--gate",
        "true" if spec["gate"] else "false",
        "--repetition",
        str(spec["repetition"]),
        "--output",
        str(record_path),
    ]
    print(
        f"[{index}/{total}] starting {spec['case']} gate={spec['gate']} "
        f"h={spec['context_depth']} r={spec['repetition']}",
        flush=True,
    )
    started = time.perf_counter()
    process = subprocess.Popen(
        command,
        cwd=PROJECT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    deadline = started + timeout
    timed_out = False
    while process.poll() is None:
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            timed_out = True
            os.killpg(process.pid, signal.SIGKILL)
            break
        try:
            process.wait(timeout=min(15.0, remaining))
        except subprocess.TimeoutExpired:
            elapsed = time.perf_counter() - started
            print(
                f"[{index}/{total}] still running {spec['case']} "
                f"gate={spec['gate']} after {elapsed:.1f}s",
                flush=True,
            )
    stdout, stderr = process.communicate()
    elapsed = time.perf_counter() - started
    stdout_path.write_text(stdout, encoding="utf-8")
    stderr_path.write_text(stderr, encoding="utf-8")

    if timed_out:
        record = failure_record(
            spec, "TIMEOUT", f"worker exceeded real process limit {timeout}s"
        )
        record["termination"] = "SIGKILL sent to fresh worker process group"
    elif process.returncode != 0:
        record = failure_record(
            spec, "PROCESS_ERROR", f"worker exited with status {process.returncode}"
        )
    elif not record_path.is_file():
        record = failure_record(spec, "INVALID_OUTPUT", "worker emitted no JSON record")
    else:
        try:
            record = json.loads(record_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            record = failure_record(spec, "INVALID_OUTPUT", str(exc))

    record.update(
        exact_command=command,
        timeout_seconds=timeout,
        process_returncode=process.returncode,
        total_process_wall_seconds=elapsed,
        stdout=stdout,
        stderr=stderr,
        stdout_log=str(stdout_path),
        stderr_log=str(stderr_path),
        record_path=str(record_path),
    )
    write_json(record_path, record)
    print(
        f"[{index}/{total}] {record['run_status']} in {elapsed:.3f}s -> {record_path}",
        flush=True,
    )
    return record


def extract_row(record: dict) -> dict:
    row = {
        "case": record["case"],
        "motifs": ",".join(record["motifs"]),
        "gate": record["gate"],
        "context_depth": record["context_depth"],
        "repetition": record["repetition"],
        "run_status": record["run_status"],
        "linear_limit": record.get("linear_limit"),
        "total_process_wall_seconds": record.get("total_process_wall_seconds"),
        "raw_record": record.get("record_path"),
    }
    analysis = record.get("analysis")
    if analysis is None:
        return row
    event_classes = analysis["composition"]["event_classes"]
    timings = analysis["diagnostics"]["seconds"]
    row.update(
        effective_configuration_sha256=analysis["configuration"]["sha256"],
        semantic_contract_sha256=analysis["contract"]["definition_sha256"],
        selected_seed_count=sum(
            origin["origin"] == "selected-seed"
            for origin in analysis["linear"]["origins"].values()
        ),
        coordinate_count=len(analysis["pjc"]["observable_coordinates"]),
        derived_event_count=len(event_classes),
        maximum_propensity_degree=max(
            event["activity_polynomial"]["degree"] for event in event_classes
        ),
        pjc_differential_order=analysis["pjc"]["operator"][
            "maximum_differential_order"
        ],
        pjc_status=analysis["pjc"]["status"],
        linear_status=analysis["linear"]["status"],
        linear_reached_dimension=len(analysis["linear"]["basis"]),
        oracle_status=analysis["word_oracle_comparison"]["status"],
        analysis_total_seconds=timings["total"],
        configuration_and_contract_seconds=timings["configuration_and_contract"],
        composition_seconds=timings["composition"],
        oracle_comparison_seconds=timings["independent_oracle_comparison"],
        linear_discovery_seconds=timings["linear_discovery"],
        pjc_seconds=timings["polynomial_jump_closure"],
    )
    return row


def summary_record(records: list[dict], protocol: dict) -> dict:
    rows = [extract_row(record) for record in records]
    completed = {}
    for row in rows:
        key = f"{row['case']}:gate={str(row['gate']).lower()}"
        completed.setdefault(key, 0)
        completed[key] += int(row["run_status"] == "COMPLETE")
    return {
        "schema_version": 1,
        "kind": "controlled-rop-scaling-summary",
        "status": (
            "COMPLETE"
            if len(rows) == protocol["requested_run_count"]
            and all(row["run_status"] == "COMPLETE" for row in rows)
            else "INCOMPLETE"
        ),
        "protocol": protocol,
        "requested_run_count": protocol["requested_run_count"],
        "recorded_run_count": len(rows),
        "actual_completed_repetitions": completed,
        "rows": rows,
    }


def report_markdown(summary: dict) -> str:
    lines = [
        "# Controlled ROP query-scaling measurements",
        "",
        "Each row is one fresh serial Python process. These are single-run diagnostic",
        "measurements, not error bars, statistical claims, or an external-tool speedup.",
        "The `abba` case includes `ABA` as well as `ABBA`; it is not the earlier",
        "dyads-plus-ABBA timing input.",
        "",
        "| Query | Gate | h | Seeds | Coordinates | Events | Degree/order | pJC | Linear (dimension) | Oracle | Process s | Config s | Composition s | Oracle s | Linear s | pJC s | Raw |",
        "|---|:---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]

    def seconds(row: dict, name: str) -> str:
        value = row.get(name)
        return "—" if value is None else f"{value:.6f}"

    for row in summary["rows"]:
        if row["run_status"] == "COMPLETE":
            linear = f"{row['linear_status']} ({row['linear_reached_dimension']})"
            order = (
                f"{row['maximum_propensity_degree']}/"
                f"{row['pjc_differential_order']}"
            )
            metrics = [
                str(row["selected_seed_count"]),
                str(row["coordinate_count"]),
                str(row["derived_event_count"]),
                order,
                row["pjc_status"],
                linear,
                row["oracle_status"],
            ]
        else:
            metrics = ["—", "—", "—", "—", row["run_status"], "—", "—"]
        raw = Path(row["raw_record"]).name if row.get("raw_record") else "—"
        lines.append(
            "| "
            + " | ".join(
                [
                    row["case"] + " (" + (row["motifs"] or "none") + ")",
                    str(row["gate"]),
                    str(row["context_depth"]),
                    *metrics,
                    seconds(row, "total_process_wall_seconds"),
                    seconds(row, "configuration_and_contract_seconds"),
                    seconds(row, "composition_seconds"),
                    seconds(row, "oracle_comparison_seconds"),
                    seconds(row, "linear_discovery_seconds"),
                    seconds(row, "pjc_seconds"),
                    f"[{raw}]({raw})" if raw != "—" else raw,
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            f"Overall status: **{summary['status']}**.",
            "",
            "`Degree/order` reports maximum propensity degree followed by maximum pJC",
            "differential order. `UNKNOWN_DIMENSION_LIMIT` is a bounded linear-worklist",
            "result and remains distinct from process failures and pJC `CLOSED`.",
            "",
        ]
    )
    return "\n".join(lines)


def run_matrix(arguments) -> None:
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    specs = schedule(
        profile=arguments.profile,
        case=arguments.case,
        gate=arguments.gate,
        repetitions=arguments.repetitions,
    )
    revision_before = git_revision()
    hashes_before = source_hashes()
    protocol = {
        "profile": arguments.profile,
        "selected_case": arguments.case,
        "gate_selection": arguments.gate,
        "requested_repetitions": arguments.repetitions,
        "requested_run_count": len(specs),
        "timeout_seconds_per_process": arguments.timeout,
        "linear_limit": LINEAR_LIMIT,
        "fresh_worker_per_configuration_and_repetition": True,
        "serial_execution": True,
        "analysis_calls_per_worker": 1,
        "precompile_before_analyze": False,
        "fixed_chemistry": "literal A/B ROP with feeds enabled",
        "event_count_expectations_role": "post-analysis sanity check only",
        "h3_difference": "includes both ABA and ABBA, unlike the earlier dyads-plus-ABBA diagnostic",
        "measurement_claim": "diagnostic observations; no error bars, speedup, or statistical claim",
    }
    run_metadata = {
        "schema_version": 1,
        "kind": "controlled-rop-scaling-run",
        "protocol": protocol,
        "schedule": specs,
        "project_revision_before": revision_before,
        "source_hashes_before": hashes_before,
        "parent_host": host_environment(),
        "parent_runtime_identity": runtime_identity(),
    }
    write_json(output / "run-metadata.json", run_metadata)
    records = [
        run_process(
            spec,
            index=index,
            total=len(specs),
            output=output,
            timeout=arguments.timeout,
        )
        for index, spec in enumerate(specs, 1)
    ]
    revision_after = git_revision()
    hashes_after = source_hashes()
    run_metadata.update(
        project_revision_after=revision_after,
        source_hashes_after=hashes_after,
        source_unchanged=hashes_before == hashes_after,
        revision_unchanged=revision_before == revision_after,
    )
    write_json(output / "run-metadata.json", run_metadata)
    summary = summary_record(records, protocol)
    summary["project_revision"] = revision_after
    summary["source_hashes"] = hashes_after
    summary["source_unchanged"] = hashes_before == hashes_after
    write_json(output / "summary.json", summary)
    fields = list(dict.fromkeys(key for row in summary["rows"] for key in row))
    with (output / "observations.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary["rows"])
    (output / "report.md").write_text(report_markdown(summary), encoding="utf-8")
    print(f"summary: {output / 'summary.json'} ({summary['status']})", flush=True)


def parse_gate(value: str) -> bool:
    return {"false": False, "true": True}[value]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    work = commands.add_parser("worker", help=argparse.SUPPRESS)
    work.add_argument("--case", choices=CASES, required=True)
    work.add_argument("--gate", type=parse_gate, choices=(False, True), required=True)
    work.add_argument("--repetition", type=int, required=True)
    work.add_argument("--output", type=Path, required=True)

    run = commands.add_parser("run")
    run.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    run.add_argument("--case", choices=CASES)
    run.add_argument("--gate", choices=("false", "true", "both"), default="both")
    run.add_argument("--repetitions", type=int, default=1)
    run.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_SECONDS)
    run.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.command == "worker":
        worker(arguments)
        return
    if arguments.repetitions < 1:
        parser.error("--repetitions must be positive")
    if arguments.timeout <= 0:
        parser.error("--timeout must be positive")
    run_matrix(arguments)


if __name__ == "__main__":
    main()
