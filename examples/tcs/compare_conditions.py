"""Compare explicit NoEdge constraints using identical overlap and valence policies.

This is an aldol-addition example runner. Composition, constraint transport,
canonicalisation and filtering are supplied by the existing implementation;
published row names are used only to label the exported results.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re

import molfoundry as mf

from mfobs.core import engine
from examples.tcs import table1


HERE = Path(__file__).resolve().parent
VARIANTS = {
    "plain": ("aldol_plain.gml", "query_unconditioned.gml"),
    "reaction_noedge": ("aldol.gml", "query_unconditioned.gml"),
    "reaction_and_query_noedge": ("aldol.gml", "query_conditioned.gml"),
}


def compare_conditions(output: Path, *, export_candidates: bool = False) -> dict:
    """Compute and export both chronological composition orders per variant."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    chemistry = json.loads((HERE / "chemistry.json").read_text(encoding="utf-8"))
    aliases = json.loads((HERE / "article-aliases.json").read_text(encoding="utf-8"))["canonical_key_sha256_to_alias"]
    names = sorted({name for pair in VARIANTS.values() for name in pair} | {"chemistry.json"})
    (output / "inputs").mkdir(exist_ok=True)
    inputs = {}
    for name in names:
        content = (HERE / name).read_bytes()
        (output / "inputs" / name).write_bytes(content)
        inputs[name] = sha256(content).hexdigest()
    report = {
        "schema_version": 1,
        "scope": {
            "overlaps": "all common overlaps, including empty; connected=False, maximum=False",
            "filter": "unchanged balanced per-atom valence and maximum C=4, H=1, O=2",
            "conditions": "explicit NoEdge transport and intrinsic edge-creation admissibility",
            "counted_terms": "signed conditional rule-isomorphism classes before observable projection",
            "excluded_scope": "no global Chem completion, C1/C2 filter, CTMC or moment hierarchy",
            "order": "first rule executes first; after minus before",
        },
        "environment": {"python": platform.python_version(), "molfoundry": version("molfoundry")},
        "input_sha256": inputs,
        "source_sha256": {
            "examples/tcs/compare_conditions.py": sha256(Path(__file__).read_bytes()).hexdigest(),
            "examples/tcs/table1.py": sha256(Path(table1.__file__).read_bytes()).hexdigest(),
            "src/mfobs/core/engine.py": sha256(Path(engine.__file__).read_bytes()).hexdigest(),
            "examples/tcs/article-aliases.json": sha256((HERE / "article-aliases.json").read_bytes()).hexdigest(),
        },
        "variants": {},
    }
    for name, (reaction_file, identity_file) in VARIANTS.items():
        reaction = mf.Rule.fromGMLString((HERE / reaction_file).read_text(encoding="utf-8"), add=False)
        identity = mf.Rule.fromGMLString((HERE / identity_file).read_text(encoding="utf-8"), add=False)
        counters, representatives, orders, contributions = {}, {}, {}, {}
        raw_exports = {}
        for order, first, second in (
            ("reaction_then_query", reaction, identity),
            ("query_then_reaction", identity, reaction),
        ):
            candidates = table1._candidates(first, second)
            decisions = [table1.chemical_rule(item.rule, chemistry) for item in candidates]
            accepted = [item for item, keep in zip(candidates, decisions) if keep]
            if export_candidates:
                # Inspection only: export the actual products before the valence
                # filter. Neither acceptance nor coefficient collection changes.
                directory = output / name / "raw"
                directory.mkdir(parents=True, exist_ok=True)
                raw_exports[order] = []
                for item, keep in zip(candidates, decisions):
                    gml = item.rule.getGMLString()
                    digest = sha256(gml.encode()).hexdigest()
                    relative = Path(name) / "raw" / (digest + ".gml")
                    (output / relative).write_text(gml, encoding="utf-8")
                    raw_exports[order].append({
                        "gml": str(relative), "accepted_by_valence_filter": keep,
                        "overlap_first_to_second": [list(pair) for pair in sorted(item.overlap)],
                    })
            counters[order] = Counter(item.rule.canonKey for item in accepted)
            contributions[order] = []
            for item in accepted:
                key = item.rule.canonKey
                representatives.setdefault(key, item.rule)
                contributions[order].append({
                    "key_sha256": sha256(key.encode()).hexdigest(),
                    "overlap_first_to_second": [list(pair) for pair in sorted(item.overlap)],
                })
            orders[order] = {
                "unfiltered": len(candidates), "filtered": len(accepted),
                "distinct_conditional_rules": len(counters[order]),
            }
        net = counters["reaction_then_query"].copy()
        net.subtract(counters["query_then_reaction"])
        rows = []
        (output / name).mkdir(exist_ok=True)
        for index, (key, rule) in enumerate(sorted(representatives.items())):
            gml = rule.getGMLString()
            # Export/reparse must preserve the conditional rule, including
            # constraints implicit in edge creation and omitted from GML.
            assert mf.Rule.fromGMLString(gml, add=False).canonKey == key
            plain_gml = re.sub(r"constrainNoEdge\s*\[[^\]]*\]", "", gml)
            plain = mf.Rule.fromGMLString(plain_gml, add=False)
            plain_hash = sha256(plain.canonKey.encode()).hexdigest()
            alias = aliases.get(plain_hash, "new-" + plain_hash[:12])
            relative = Path(name) / f"composite-{index:02}.gml"
            (output / relative).write_text(gml, encoding="utf-8")
            rows.append({
                "key_sha256": sha256(key.encode()).hexdigest(),
                "plain_class": alias,
                "after": counters["reaction_then_query"][key],
                "before": counters["query_then_reaction"][key],
                "net_coefficient": net[key], "gml": str(relative),
                "constraints": [{'source': c.source, 'target': c.target} for c in rule.matchConstraintsNoEdge],
            })
        report["variants"][name] = {
            "inputs": {"reaction": reaction_file, "identity": identity_file},
            "orders": orders, "contributions": contributions,
            "distinct_rules_union": len(representatives),
            "nonzero_signed_terms": sum(value != 0 for value in net.values()),
            "signed_term_absolute_multiplicity": sum(abs(value) for value in net.values()),
            "rows": rows,
        }

        if export_candidates:
            report["variants"][name]["raw_candidates"] = raw_exports

    base = report["variants"]["reaction_noedge"]
    conditioned = report["variants"]["reaction_and_query_noedge"]
    removed = {}
    for order in base["orders"]:
        overlap_key = lambda row: tuple(tuple(pair) for pair in row["overlap_first_to_second"])
        original = {overlap_key(row): row for row in base["contributions"][order]}
        kept = {overlap_key(row): row for row in conditioned["contributions"][order]}
        assert kept.keys() <= original.keys(), "Adding NoEdge must only restrict overlaps"
        row_by_key = {row["key_sha256"]: row for row in base["rows"]}
        removed[order] = [
            {**original[key], "plain_class": row_by_key[original[key]["key_sha256"]]["plain_class"]}
            for key in sorted(original.keys() - kept.keys())
        ]
    report["comparison"] = {
        "removed_contributions": removed,
        "reaction_constraint_held_fixed": True,
        "selected_r180": {
            name: [{field: row[field] for field in ("gml", "after", "before", "net_coefficient", "constraints")}
                   for row in data["rows"] if row["plain_class"] == "r180"]
            for name, data in report["variants"].items()
        },
    }
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    report = compare_conditions(args.out_dir)
    for name, result in report["variants"].items():
        print(name, result["orders"], "nonzero_signed_terms", result["nonzero_signed_terms"])


if __name__ == "__main__":
    main()
