#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from mfobs.core.source_identity import FRAGGER_SNAPSHOT, require_source_identity

FRAGGER_COMMIT = "0d9595c5d6e8918fc605f7d65778568a74f5e2d4"
JAR_HASHES = {
    "scala-compiler-2.13.2.jar": "b04da83111decff5b36e87c74d73f5e2231b63ce605f2602e862bae917d7761a",
    "scala-library-2.13.2.jar": "6a1c9846feb8f5d7d33b91e61e056d272215a79747a2b7b5e5ece690bf151d5d",
    "scala-reflect-2.13.2.jar": "21a7288874f186d6049e8b14a986011ffd96cb587f206eccf946624f2d476815",
}
LENGTHS = (2,4,6,8,10,12,14,16)
FAMILIES = ("alt","same")

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(command, *, cwd=None, timeout=60, check=True):
    return subprocess.run(
        [str(x) for x in command],
        cwd=cwd,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=check,
    )

def parse_record(stdout: str, prefix: str) -> dict:
    rows = [line[len(prefix):] for line in stdout.splitlines() if line.startswith(prefix)]
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one {prefix} record, found {len(rows)}")
    return json.loads(rows[0])

def verify_inputs(args):
    checkout = args.fragger_checkout.resolve()
    wheel = args.wheel.resolve() if args.wheel is not None else None
    java_version = run([args.java,"-version"],check=False).stderr.strip()
    if '"11.' not in java_version and 'version 11' not in java_version:
        raise RuntimeError("publication benchmark requires a Java 11 runtime for the pinned Scala 2.13.2 core")
    fragger_identity = require_source_identity(
        checkout,
        expected_commit=FRAGGER_COMMIT,
        snapshot_name=FRAGGER_SNAPSHOT,
        git_paths=("src/main",),
    )
    compiler = checkout / "compiler"
    observed = {
        name: sha256(compiler / name)
        for name in JAR_HASHES
    }
    if observed != JAR_HASHES:
        raise RuntimeError(f"Scala compiler jar hashes differ: {observed}")
    return checkout, wheel, compiler, fragger_identity

def compile_fragger(checkout: Path, compiler: Path, java: Path, out: Path):
    classes = out / "fragger-classes"
    classes.mkdir()
    sources = sorted((checkout / "src/main/scala/graph-rewriting").glob("*.scala"))
    harness = HERE / "ABScalingHarness.scala"
    command = [
        java, "-Xmx2g", "-cp", str(compiler / "*"),
        "scala.tools.nsc.Main", "-usejavacp", "-d", classes,
        *sources, harness,
    ]
    start = time.perf_counter()
    result = run(command,cwd=checkout,timeout=180)
    elapsed = time.perf_counter() - start
    (out / "fragger-compile.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    return classes, command, elapsed

def invoke_molfoundry(args, wheel: Path | None, family: str, n: int, out: Path) -> dict:
    command = [
        args.python, HERE / "molfoundry_benchmark.py",
        family, str(n),
        "--warmups", str(args.warmups),
        "--repetitions", str(args.repetitions),
        "--molfoundry-version", args.molfoundry_version,
    ]
    if wheel is not None:
        command += ["--wheel", wheel]
    start = time.perf_counter()
    result = run(command,cwd=ROOT,timeout=args.timeout)
    process_seconds = time.perf_counter() - start
    log = out / f"molfoundry-{family}-{n}.log"
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    meta = parse_record(result.stdout,"META=")
    record = parse_record(result.stdout,"RESULT=")
    record["runtime_meta"] = meta
    record["process_seconds"] = process_seconds
    record["log"] = log.relative_to(out).as_posix()
    return record

def invoke_fragger(args, checkout: Path, compiler: Path, classes: Path,
                   family: str, n: int, out: Path) -> dict:
    classpath = str(classes) + ":" + str(compiler / "*")
    command = [
        args.java, "-Xmx4g", "-cp", classpath,
        "hz.ricardo.graph_rewriting.ABScalingHarness",
        family, str(n), str(args.warmups), str(args.repetitions),
    ]
    start = time.perf_counter()
    result = run(command,cwd=checkout,timeout=args.timeout)
    process_seconds = time.perf_counter() - start
    log = out / f"fragger-{family}-{n}.log"
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    record = parse_record(result.stdout,"RESULT=")
    record["process_seconds"] = process_seconds
    record["log"] = log.relative_to(out).as_posix()
    return record

def environment(args, checkout: Path, wheel: Path | None, compiler: Path,
                compile_command, compile_seconds, molfoundry_meta: dict,
                fragger_identity: dict):
    java_version = run([args.java,"-version"],check=False).stderr.strip()
    python_version = run([args.python,"--version"]).stdout.strip()
    lscpu = run(["lscpu"],check=False).stdout if shutil.which("lscpu") else None
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": python_version,
        "java": java_version,
        "lscpu": lscpu,
        "fragger_commit": FRAGGER_COMMIT,
        "fragger_checkout": str(checkout),
        "fragger_source_identity": fragger_identity,
        "scala_jar_sha256": {name: sha256(compiler/name) for name in JAR_HASHES},
        "molfoundry_wheel": str(wheel) if wheel is not None else None,
        "molfoundry_wheel_sha256": sha256(wheel) if wheel is not None else None,
        "molfoundry_runtime": molfoundry_meta,
        "benchmark_source_sha256": {
            path.name: sha256(path)
            for path in (
                HERE/"ABScalingHarness.scala",
                HERE/"molfoundry_benchmark.py",
                HERE/"run_benchmark.py",
            )
        },
        "fragger_compile_command": [str(x) for x in compile_command],
        "fragger_compile_seconds": compile_seconds,
    }

def summarize(records, families, lengths):
    by_key = {}
    for row in records:
        by_key.setdefault((row["family"],row["n"]),{})[row["tool"]] = row
    rows = []
    for family in families:
        for n in lengths:
            pair = by_key.get((family,n),{})
            if set(pair) != {"molfoundry","fragger-core"}:
                continue
            mf = pair["molfoundry"]
            fr = pair["fragger-core"]
            rows.append({
                "family": family,
                "n": n,
                "word": mf["word"],
                "molfoundry_median_seconds": mf["median_seconds"],
                "fragger_median_seconds": fr["median_seconds"],
                "fragger_over_molfoundry": fr["median_seconds"]/mf["median_seconds"],
                "molfoundry_raw_composites": mf["raw_composites"],
                "fragger_raw_composites": fr["raw_composites"],
            })
    return rows

def verify_counts(records):
    expected = {
        (f"append_{previous}_{symbol}", order)
        for symbol in ("A", "B") for previous in ("R", "A", "B")
        for order in ("observable-left", "observable-right")
    }
    counts = {}
    for record in records:
        rows = record["per_order"]
        per_order = {(row["rule"], row["order"]): row["count"] for row in rows}
        if len(rows) != 12 or set(per_order) != expected:
            raise ValueError(f"{record['tool']}: expected twelve distinct ordered compositions")
        if sum(per_order.values()) != record["raw_composites"]:
            raise ValueError(f"{record['tool']}: total differs from per-order counts")
        counts[record["tool"]] = per_order
    if set(counts) != {"molfoundry", "fragger-core"}:
        raise ValueError("missing tool result")
    mismatches = [
        f"{rule}/{order}: MolFoundry={counts['molfoundry'][rule, order]}, "
        f"FRAGGER={counts['fragger-core'][rule, order]}"
        for rule, order in sorted(expected)
        if counts["molfoundry"][rule, order] != counts["fragger-core"][rule, order]
    ]
    if mismatches:
        raise ValueError("; ".join(mismatches))


def write_csv(rows, path: Path):
    columns = [
        "family","n","word","molfoundry_median_seconds","fragger_median_seconds",
        "fragger_over_molfoundry","molfoundry_raw_composites","fragger_raw_composites",
    ]
    lines = [",".join(columns)]
    for row in rows:
        fields = []
        for col in columns:
            value = str(row[col])
            if "," in value or '"' in value:
                value = '"' + value.replace('"','""') + '"'
            fields.append(value)
        lines.append(",".join(fields))
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")

def write_checkpoint(out: Path, args, checkout: Path, wheel: Path | None,
                     compiler: Path, compile_command, compile_seconds,
                     fragger_identity: dict, records: list[dict],
                     status: str, failure: dict | None = None):
    rows = summarize(records,args.families,args.lengths)
    mf_rows = [row for row in records if row.get("tool") == "molfoundry"]
    manifest = {
        "schema_version": 1,
        "kind": "ab-composition-scaling-benchmark",
        "status": status,
        "protocol": {
            "families": args.families,
            "requested_molfoundry_version": args.molfoundry_version,
            "lengths": args.lengths,
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "timeout_seconds_per_process": args.timeout,
            "serial": True,
            "alternating_tool_order": True,
            "primary_timing_scope": "in-process composition enumeration after warm-up; excludes process startup/import/compilation",
            "raw_count_scope": "returned composites from 12 ordered rule/observable compositions, before canonical deduplication/cancellation",
        },
        "environment": (
            environment(
                args,checkout,wheel,compiler,compile_command,compile_seconds,
                mf_rows[0]["runtime_meta"],fragger_identity,
            )
            if mf_rows else None
        ),
        "failure": failure,
        "observations": records,
        "comparison": rows,
    }
    (out/"results.json").write_text(json.dumps(manifest,indent=2)+"\n", encoding="utf-8")
    write_csv(rows,out/"results.csv")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragger-checkout", type=Path, required=True)
    parser.add_argument("--wheel", type=Path)
    parser.add_argument("--molfoundry-version", default="0.1.9")
    parser.add_argument("--java", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=list(FAMILIES))
    parser.add_argument("--lengths", nargs="+", type=int, default=list(LENGTHS))
    args = parser.parse_args()

    if args.repetitions < 3:
        parser.error("publication run requires at least 3 timed repetitions")
    if any(n < 2 for n in args.lengths):
        parser.error("all lengths must be >=2")

    out = args.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    checkout,wheel,compiler,fragger_identity = verify_inputs(args)
    classes,compile_command,compile_seconds = compile_fragger(
        checkout,compiler,args.java.resolve(),out
    )

    records = []
    case_index = 0
    write_checkpoint(
        out,args,checkout,wheel,compiler,compile_command,compile_seconds,
        fragger_identity,records,"RUNNING",
    )
    for family in args.families:
        for n in args.lengths:
            case_records = []
            # Alternate tool order to reduce monotone thermal/order bias.
            order = ("molfoundry","fragger") if case_index % 2 == 0 else ("fragger","molfoundry")
            for tool in order:
                try:
                    if tool == "molfoundry":
                        row = invoke_molfoundry(args,wheel,family,n,out)
                    else:
                        row = invoke_fragger(args,checkout,compiler,classes,family,n,out)
                except subprocess.TimeoutExpired as exc:
                    failure = {
                        "kind": "TIMEOUT",
                        "tool": tool,
                        "family": family,
                        "n": n,
                        "timeout_seconds": args.timeout,
                        "command": [str(x) for x in exc.cmd],
                    }
                    write_checkpoint(
                        out,args,checkout,wheel,compiler,compile_command,compile_seconds,
                        fragger_identity,records,"TIMEOUT",failure,
                    )
                    print(
                        "TIMEOUT",tool,family,n,f"after {args.timeout}s",
                        flush=True,
                    )
                    return
                records.append(row)
                case_records.append(row)
                write_checkpoint(
                    out,args,checkout,wheel,compiler,compile_command,compile_seconds,
                    fragger_identity,records,"RUNNING",
                )
                print(
                    tool, family, n,
                    f'{row["median_seconds"]:.6g}s',
                    row["raw_composites"],
                    flush=True,
                )
            try:
                verify_counts(case_records)
            except ValueError as exc:
                failure = {"kind": "COUNT_MISMATCH", "family": family, "n": n,
                           "message": str(exc)}
                write_checkpoint(
                    out,args,checkout,wheel,compiler,compile_command,compile_seconds,
                    fragger_identity,records,"FAIL",failure,
                )
                raise SystemExit(f"COUNT_MISMATCH {family} {n}: {exc}") from exc
            case_index += 1

    write_checkpoint(
        out,args,checkout,wheel,compiler,compile_command,compile_seconds,
        fragger_identity,records,"PASS",
    )
    print("PASS",out/"results.json")

if __name__ == "__main__":
    main()
