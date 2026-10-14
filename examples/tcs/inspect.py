"""Inspect the local TCS replay and its computed composition table.

The published table remains a comparison oracle.  The displayed table and its
per-row GML are recomputed from the visible rules and chemistry declaration.
A commutator is exported as GML terms plus signed coefficients in JSON because
a linear combination of rules is not itself one GML rule.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from .table1 import export_computed_table, render_computed_table


HERE = Path(__file__).resolve().parent
ARTICLE_TABLE = HERE / "article-table.json"


def load_article_table() -> dict:
    table = json.loads(ARTICLE_TABLE.read_text(encoding="utf-8"))
    rows = table["rows"]
    if len(rows) != 24:
        raise ValueError("TCS Table 1 transcription must contain 24 rows")
    first = sum(row["r_plus_after_r_id"] for row in rows)
    second = sum(row["r_id_after_r_plus"] for row in rows)
    if (first, second) != (12, 18):
        raise ValueError("TCS Table 1 transcription must contain 12/18 contributions")
    if table["manual_only_rules"] != ["rA", "rB", "rC"]:
        raise ValueError("unexpected manual-only rule list")
    return table


def render_published_table(table: Mapping | None = None) -> str:
    """Render the publication oracle; this is never the CLI's computed table."""
    table = dict(table or load_article_table())
    headings = ("rule", "|r+ * rid|", "|rid * r+|", "C1", "C2")
    body = [
        (
            row["rule"],
            str(row["r_plus_after_r_id"] or ""),
            str(row["r_id_after_r_plus"] or ""),
            "X" if row["violates_C1"] else "",
            "X" if row["violates_C2"] else "",
        )
        for row in table["rows"]
    ]
    widths = [max(len(headings[i]), *(len(row[i]) for row in body)) for i in range(5)]
    line = "  ".join(headings[i].ljust(widths[i]) for i in range(5))
    separator = "  ".join("-" * width for width in widths)
    rows = ["  ".join(row[i].ljust(widths[i]) for i in range(5)) for row in body]
    source = table["source"]
    return "\n".join([
        "Published TCS composition counts - comparison oracle, not replay output",
        f"Source: {source['arxiv']}, {source['location']}",
        line,
        separator,
        *rows,
        "Totals: 12 and 18 contributions; rA, rB, rC were manual-only in the article.",
    ])


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _host_gml(host: Mapping) -> str:
    if host.get("gml"):
        return host["gml"]
    lines = ["graph ["]
    lines.extend(
        f'  node [ id {node["id"]} label {json.dumps(node["label"])} ]'
        for node in host["nodes"]
    )
    lines.extend(
        f'  edge [ source {edge["source"]} target {edge["target"]} '
        f'label {json.dumps(edge["label"])} ]'
        for edge in host["edges"]
    )
    return "\n".join([*lines, "]", ""])


def _prepare_output(output: Path, *, force: bool) -> None:
    if output.exists() and any(output.iterdir()) and not force:
        raise FileExistsError(f"{output} is not empty; pass --force to replace named files")
    output.mkdir(parents=True, exist_ok=True)


def _export_selected(record: Mapping, output: Path) -> dict:
    inputs_dir = output / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    for name, gml in record["inputs"].items():
        (inputs_dir / f"{name}.gml").write_text(gml, encoding="utf-8")

    witness_file = output / "host-witness.gml"
    witness_file.write_text(_host_gml(record["host_witness"]), encoding="utf-8")

    selected_manifest = {}
    commutator_manifest = {}
    for variant, data in record["variants"].items():
        selected_manifest[variant] = []
        commutator_manifest[variant] = []
        for side in ("before", "after"):
            composite = data["selected_composites"][side]
            relative = Path("selected-composites") / variant / f"{side}.gml"
            path = output / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(composite["gml"], encoding="utf-8")
            selected_manifest[variant].append({
                "order": side,
                "gml": str(relative),
                "canonical_key_sha256": composite["key_sha256"],
                "overlap_first_to_second": composite["overlap_first_to_second"],
            })
            coefficient = data["selected_net_commutator_coefficients"][side]
            if coefficient:
                term_relative = Path("commutators") / variant / f"{side}.gml"
                term = output / term_relative
                term.parent.mkdir(parents=True, exist_ok=True)
                term.write_text(composite["gml"], encoding="utf-8")
                commutator_manifest[variant].append({
                    "coefficient": coefficient,
                    "gml": str(term_relative),
                    "canonical_key_sha256": composite["key_sha256"],
                })
    _write_json(output / "selected-composites.json", selected_manifest)
    _write_json(output / "selected-commutators.json", {
        "meaning": "Each list is a signed linear combination of the referenced GML rules.",
        "terms": commutator_manifest,
    })
    return {
        "inputs": len(record["inputs"]),
        "selected_composites": sum(map(len, selected_manifest.values())),
        "selected_commutator_terms": sum(map(len, commutator_manifest.values())),
        "host_witness": witness_file.relative_to(output).as_posix(),
    }


def _group_rules(rules: Sequence) -> tuple[Counter, dict]:
    counts = Counter(rule.canonKey for rule in rules)
    representatives = {}
    for rule in rules:
        representatives.setdefault(rule.canonKey, rule)
    return counts, representatives


def _export_all(products: Mapping, output: Path) -> dict:
    summary = {}
    for variant, orders in products.items():
        grouped = {side: _group_rules(orders[side]) for side in ("before", "after")}
        summary[variant] = {"compositions": {}, "commutator_terms": 0}
        for side in ("before", "after"):
            counts, representatives = grouped[side]
            rows = []
            for index, key in enumerate(sorted(counts), 1):
                relative = Path("all-compositions") / variant / side / f"term-{index:04d}.gml"
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(representatives[key].getGMLString(), encoding="utf-8")
                rows.append({
                    "multiplicity": counts[key], "gml": str(relative),
                    "canonical_key_sha256": sha256(key.encode()).hexdigest(),
                })
            _write_json(output / "all-compositions" / variant / f"{side}.json", rows)
            summary[variant]["compositions"][side] = len(rows)

        before_counts, before_reps = grouped["before"]
        after_counts, after_reps = grouped["after"]
        rows = []
        for key in sorted(set(before_counts) | set(after_counts)):
            coefficient = after_counts[key] - before_counts[key]
            if not coefficient:
                continue
            representative = after_reps.get(key) or before_reps[key]
            index = len(rows) + 1
            relative = Path("all-commutators") / variant / f"term-{index:04d}.gml"
            path = output / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(representative.getGMLString(), encoding="utf-8")
            rows.append({
                "coefficient": coefficient, "gml": str(relative),
                "canonical_key_sha256": sha256(key.encode()).hexdigest(),
            })
        _write_json(output / "all-commutators" / variant / "manifest.json", {
            "meaning": "Signed commutator after minus before; each term references one GML rule.",
            "terms": rows,
        })
        summary[variant]["commutator_terms"] = len(rows)
    return summary


def export_gml(record: Mapping, output: Path, *, products: Mapping | None = None,
               force: bool = False) -> dict:
    """Export transparent GML files and signed manifests for a TCS replay."""
    _prepare_output(output, force=force)
    table = load_article_table()
    _write_json(output / "published-table-oracle.json", table)
    computed_table = export_computed_table(output)
    summary = {
        "scope": (
            "Local labelled-graph replay plus a fresh computation of the plain aldol "
            "composition table from declared chemistry constraints. The published table "
            "is retained only as a comparison oracle."
        ),
        "computed_table": computed_table,
        "selected": _export_selected(record, output),
        "all_computed": _export_all(products, output) if products is not None else None,
    }
    _write_json(output / "manifest.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path, help="JSON produced by the TCS replay")
    parser.add_argument(
        "--gml-dir", type=Path,
        help="Export the computed table row GML, selected terms, and manifests",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-table", action="store_true")
    args = parser.parse_args()
    record = json.loads(args.record.read_text(encoding="utf-8"))
    if not args.no_table:
        print(render_computed_table())
    if args.gml_dir:
        summary = export_gml(record, args.gml_dir, force=args.force)
        print(f"\nGML export: {args.gml_dir}")
        print(json.dumps(summary["selected"], sort_keys=True))


if __name__ == "__main__":
    main()
