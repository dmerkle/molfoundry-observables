#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import statistics
import sys
import time

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import molfoundry
from molfoundry import Rule
from mfobs.core.typed_incidence import grouped_typed_incidence_compositions

SPECIFICATION = {
    "name": "ab-directed-multigraph-incidence-v1",
    "objectLabels": ["R", "A", "B", "T"],
    "edgeNodeLabels": ["Next"],
    "sourceBond": "src",
    "targetBond": "tgt",
    "atMostOneOutgoing": False,
    "atMostOneIncomingTypes": [],
    "forbidSuccessorAfterTypes": [],
    "simpleDirected": False,
    "simpleUnderlying": False,
    "acyclicDirected": False,
    "edgeEndpointTypes": {
        "Next": {
            "source": ["R", "A", "B", "T"],
            "target": ["R", "A", "B", "T"],
        }
    },
    "relationBonds": {},
    "relationIncidentBounds": [],
}

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def incidence_path(word: str) -> str:
    nodes = [f'node [ id {i} label "{symbol}" ]' for i, symbol in enumerate(word)]
    edges = []
    next_id = len(word)
    for i in range(len(word) - 1):
        edge_id = next_id
        next_id += 1
        nodes.append(f'node [ id {edge_id} label "Next" ]')
        edges.append(f'edge [ source {i} target {edge_id} label "src" ]')
        edges.append(f'edge [ source {edge_id} target {i + 1} label "tgt" ]')
    return "\n".join(nodes + edges)

def identity(word: str) -> Rule:
    return Rule.fromGMLString(
        "rule [\n context [\n" + incidence_path(word) + "\n ]\n]",
        name="O_" + word,
        add=False,
    )

def append(previous: str, symbol: str) -> Rule:
    return Rule.fromGMLString(
        f"""rule [
 context [
  node [ id 0 label "{previous}" ]
  node [ id 1 label "T" ]
 ]
 left [
  node [ id 2 label "Next" ]
  edge [ source 0 target 2 label "src" ]
  edge [ source 2 target 1 label "tgt" ]
 ]
 right [
  node [ id 3 label "{symbol}" ]
  node [ id 4 label "Next" ]
  node [ id 5 label "Next" ]
  edge [ source 0 target 4 label "src" ]
  edge [ source 4 target 3 label "tgt" ]
  edge [ source 3 target 5 label "src" ]
  edge [ source 5 target 1 label "tgt" ]
 ]
]""",
        name=f"append_{previous}_{symbol}",
        add=False,
    )

RULES = [
    append(previous, symbol)
    for symbol in ("A", "B")
    for previous in ("R", "A", "B")
]

def expressions(observable: Rule):
    rows = []
    for rule in RULES:
        rows.append((rule.name, "observable-left", observable, rule))
        rows.append((rule.name, "observable-right", rule, observable))
    return rows

def compose_count(first: Rule, second: Rule) -> tuple[int, dict]:
    diagnostics, counts, _ = grouped_typed_incidence_compositions(
        first,
        second,
        specification=SPECIFICATION,
        connected=False,
        maximum=False,
        include_empty=True,
    )
    count = sum(counts.values())
    if count != diagnostics.accepted_raw:
        raise RuntimeError("typed-incidence raw multiplicity mismatch")
    return count, diagnostics.record()

def enumerate_raw(rows):
    per_order = []
    total = 0
    for rule, order, first, second in rows:
        count, diagnostics = compose_count(first, second)
        per_order.append({
            "rule": rule,
            "order": order,
            "count": count,
            "diagnostics": diagnostics,
        })
        total += count
    return {"raw_composites": total, "per_order": per_order}

def word_for(family: str, n: int) -> str:
    if family == "same":
        return "A" * n
    return ("AB" * ((n + 1) // 2))[:n]

def benchmark(family: str, n: int, warmups: int, repetitions: int) -> dict:
    word = word_for(family, n)
    rows = expressions(identity(word))

    reference = None
    for _ in range(warmups):
        reference = enumerate_raw(rows)

    times = []
    result = None
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        result = enumerate_raw(rows)
        stop = time.perf_counter_ns()
        if reference is not None and result != reference:
            raise RuntimeError("non-deterministic raw composition result")
        reference = result
        times.append((stop - start) / 1e9)

    return {
        "tool": "molfoundry",
        "representation": "directed-multigraph-via-typed-incidence",
        "family": family,
        "word": word,
        "n": n,
        **result,
        "times_seconds": times,
        "median_seconds": statistics.median(times),
        "min_seconds": min(times),
        "max_seconds": max(times),
    }

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("family", choices=("alt", "same"))
    parser.add_argument("n", type=int)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--wheel", type=Path)
    parser.add_argument("--molfoundry-version", default="0.1.9")
    args = parser.parse_args()

    if args.n < 2:
        parser.error("n must be at least 2")
    if args.warmups < 0 or args.repetitions < 1:
        parser.error("invalid warmup/repetition count")

    installed_version = importlib.metadata.version("molfoundry")
    if installed_version != args.molfoundry_version:
        raise SystemExit(f"expected molfoundry {args.molfoundry_version}, found {installed_version}")

    package_root = Path(molfoundry.__file__).resolve().parent
    engine = package_root / "bin" / "molfoundry.js"
    engine_sha = sha256(engine)
    wheel_sha = None
    if args.wheel is not None:
        wheel_sha = sha256(args.wheel.resolve())

    global_rows = expressions(identity("ABABA"))
    for _ in range(4):
        enumerate_raw(global_rows)

    selected_files = {
        name: sha256(package_root / name)
        for name in ("_engine.py", "rule.py", "rulecomposition.py", "bin/molfoundry.js")
        if (package_root / name).is_file()
    }
    typed_incidence = PROJECT_ROOT / "src/mfobs/core/typed_incidence.py"
    composition = PROJECT_ROOT / "src/mfobs/core/composition.py"
    print("META=" + json.dumps({
        "molfoundry_version": importlib.metadata.version("molfoundry"),
        "package_root": str(package_root),
        "selected_file_sha256": selected_files,
        "engine_sha256": engine_sha,
        "wheel_sha256": wheel_sha,
        "typed_incidence_sha256": sha256(typed_incidence),
        "composition_sha256": sha256(composition),
        "python": sys.version,
    }, sort_keys=True))
    print("RESULT=" + json.dumps(
        benchmark(args.family, args.n, args.warmups, args.repetitions),
        sort_keys=True,
    ))

if __name__ == "__main__":
    main()
