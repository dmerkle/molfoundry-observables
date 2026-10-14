#!/usr/bin/env python3
"""Reproduce declarative-family closures, negative run, and independent checks."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="Skip the degree-10 and degree-11 scaling cases")
    parser.add_argument("--summarize-only", action="store_true", help="Read existing result/check files; do not rerun synthesis")
    parser.add_argument("--out", type=Path, default=ROOT / "out/gated-path-checks",
                        help="Output root; relative paths are relative to the project root")
    args = parser.parse_args()
    output_root = (args.out if args.out.is_absolute() else ROOT / args.out).resolve()
    cases = [("model.json", "gated-path", []),
             ("model-ungated.json", "ungated-path", ["--ungated"]),
             ("model-degree2.json", "degree2", []),
             ("model-degree4.json", "degree4", ["--max-n", "5"])]
    if not args.quick:
        cases.append(("model-degree10.json", "degree10", ["--max-n", "3"]))
        cases.append(("model-degree11.json", "degree11", ["--max-n", "3"]))
    cases.append(("model-feedback.json", "feedback-path", ["--feedback", "--partial", "--max-n", "5"]))
    expected_dimensions = {"gated-path": 3, "ungated-path": 2, "degree2": 2,
                           "degree4": 6, "degree10": 18, "degree11": 20}
    summary = []
    for filename, output, check_options in cases:
        model = ROOT / "examples/gated-path" / filename
        destination = output_root / output
        compressed_destination = output_root / f"{output}-compressed"
        command = [sys.executable, "run.py", "run", str(model), "--out", str(destination), "--progress", "--force"]
        if output == "feedback-path":
            command += ["--max-expansions", "5"]
        if not args.summarize_only:
            completed = subprocess.run(command, cwd=ROOT, check=False)
            assert completed.returncode == (2 if output == "feedback-path" else 0), command
            subprocess.run([sys.executable, "scripts/validation/verify_gated_path.py", str(model),
                            str(destination / "result.json"), "--out", str(destination / "independent-checks.json"),
                            *check_options], cwd=ROOT, check=True)
        result = json.loads((destination / "result.json").read_text(encoding="utf-8"))
        checks = json.loads((destination / "independent-checks.json").read_text(encoding="utf-8"))
        assert checks["passed"] and checks["basis_size"] == len(result["basis"])
        compressed_dimension = None
        if result["closed"]:
            if not args.summarize_only:
                subprocess.run([sys.executable, "run.py", "compress", str(destination / "result.json"),
                                "--out", str(compressed_destination), "--force"], cwd=ROOT, check=True)
            compressed_path = compressed_destination / "result.json"
            if not compressed_path.is_file():
                raise SystemExit(f"Missing compressed result {compressed_path}; summarize-only never invents or recomputes it")
            compressed = json.loads(compressed_path.read_text(encoding="utf-8"))
            source_snapshot = json.dumps(result, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
            assert compressed["source_result_sha256"] == hashlib.sha256(source_snapshot.encode()).hexdigest(), "Compressed result refers to a different source run"
            assert compressed["closed"] and compressed["dimension"] == expected_dimensions[output]
            compressed_dimension = compressed["dimension"]
        summary.append({"case": output, "status": result["status"], "basis_size": len(result["basis"]),
                        "compressed_dimension": compressed_dimension,
                        "rounds": [r["after"] for r in result.get("rounds", [])],
                        "raw_compositions": result["raw_compositions"], "synthesis_seconds": result["seconds"],
                        "symbolic_generator_identities_all_n": checks["symbolic_generator_identities_all_n"],
                        "finite_graph_generator_identities": checks["independent_graph_generator_identities"],
                        "specialized_krylov_rank": checks.get("independent_minimality_diagnostic", {}).get("rank")})
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
