"""Recompute the aldol composition table and compare it with TCS Table 1.

Nothing in the published table is used to construct a composite rule.  The
rules, overlap policy, and chemical constraints are executable inputs.  The
published table and historical rule aliases are comparison or presentation
oracles only.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

import molfoundry as mf

from mfobs.core.composition import witnessed_compositions
from mfobs.core.graph import edge_rows
from mfobs.core.gml import rule_side
from examples.tcs.inputs import chemical_rule, _side


HERE = Path(__file__).resolve().parent
CHEMISTRY = HERE / "chemistry.json"
PUBLISHED = HERE / "article-table.json"
ALIASES = HERE / "article-aliases.json"


@dataclass(frozen=True)
class Candidate:
    rule: object
    overlap: tuple[tuple[int, int], ...]
    edge_reflecting: bool


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))



def convenience_flags(rule, chemistry: Mapping) -> dict[str, bool]:
    """Evaluate declared local star predicates on both rule interfaces."""
    sides = [
        _side(rule, side, chemistry["bond_orders"])
        for side in ("left", "right")
    ]

    def violates(specification: Mapping) -> bool:
        kind = specification["kind"]
        center = specification["center_label"]
        minimum = int(specification["minimum"])
        for labels, incident, _ in sides:
            for node, label in labels.items():
                if label != center:
                    continue
                if kind == "incident_edge_count":
                    count = sum(
                        edge == specification["edge_label"]
                        for edge, _ in incident[node]
                    )
                elif kind == "neighbour_label_count":
                    count = sum(
                        neighbour == specification["neighbour_label"]
                        for _, neighbour in incident[node]
                    )
                else:
                    raise ValueError(f"unsupported reported constraint kind: {kind}")
                if count >= minimum:
                    return True
        return False

    return {
        f"violates_{name}": violates(specification)
        for name, specification in chemistry["reported_constraints"].items()
    }


def _candidates(first, second) -> list[Candidate]:
    right, left = rule_side(first, 'right'), rule_side(second, 'left')
    right_edges = {(u,v): label for u,v,label in edge_rows(right)}
    left_edges = {(u,v): label for u,v,label in edge_rows(left)}
    def reflects(overlap):
        from itertools import combinations
        return all(right_edges.get(tuple(sorted((a,c)))) ==
                   left_edges.get(tuple(sorted((b,d))))
                   for (a,b),(c,d) in combinations(overlap,2))
    return [Candidate(rule, overlap, reflects(overlap))
            for rule, overlap in witnessed_compositions(first, second)]


def _alias_map() -> dict[str, str]:
    return _load_json(ALIASES)["canonical_key_sha256_to_alias"]


def normalized_published_rows(publication=None) -> list[dict]:
    """Merge verified plain-rule aliases; preserve the literal transcription."""
    audit = _load_json(HERE / "article-manual-aliases.json")
    aliases = _alias_map()
    for row in audit["exact_overlap_evidence"]:
        digest = sha256(row["canonical_key"].encode()).hexdigest()
        if aliases.get(digest) != audit["mapping"][row["article_name"]]:
            raise ValueError("Manual alias evidence disagrees with historical class identity")
    merged = {}
    for row in (publication or _load_json(PUBLISHED))["rows"]:
        name = audit["mapping"].get(row["rule"], row["rule"])
        if name not in merged:
            merged[name] = {**row, "rule": name, "published_names": [row["rule"]]}
        else:
            target = merged[name]
            target["published_names"].append(row["rule"])
            for field in ("r_plus_after_r_id", "r_id_after_r_plus"):
                target[field] += row[field]
            # Do not infer agreement for contradictory published annotations.
            for field in ("violates_C1", "violates_C2"):
                if target[field] != row[field]:
                    raise ValueError("Aliased published rows have conflicting constraint annotations")
    return list(merged.values())


def _published_rows() -> tuple[list[str], dict[str, dict]]:
    rows = normalized_published_rows()
    return [row["rule"] for row in rows], {row["rule"]: row for row in rows}


def _name(key: str, aliases: Mapping[str, str]) -> str:
    digest = sha256(key.encode()).hexdigest()
    return aliases.get(digest, f"new-{digest[:12]}")


def _analyze() -> tuple[dict, dict[str, object]]:
    chemistry = _load_json(CHEMISTRY)
    aliases = _alias_map()
    published_order, published = _published_rows()
    query = mf.Rule.fromGMLString((HERE / "query_unconditioned.gml").read_text(encoding="utf-8"), add=False)
    reaction = mf.Rule.fromGMLString((HERE / "aldol_plain.gml").read_text(encoding="utf-8"), add=False)

    # Multiplication is written in rule-algebra order.  The first tuple below
    # therefore corresponds to the first numerical column of Table 1.
    orders = {
        "r_plus_after_r_id": (reaction, query),
        "r_id_after_r_plus": (query, reaction),
    }
    counters: dict[str, Counter] = {}
    automatic: dict[str, Counter] = {}
    representatives: dict[str, object] = {}
    non_reflecting: dict[str, list[dict]] = {}

    for order, (first, second) in orders.items():
        accepted = [
            candidate for candidate in _candidates(first, second)
            if chemical_rule(candidate.rule, chemistry)
        ]
        counters[order] = Counter(item.rule.canonKey for item in accepted)
        automatic[order] = Counter(
            item.rule.canonKey for item in accepted if item.edge_reflecting
        )
        non_reflecting[order] = []
        for item in accepted:
            key = item.rule.canonKey
            representatives.setdefault(key, item.rule)
            if not item.edge_reflecting:
                non_reflecting[order].append({
                    "class": _name(key, aliases),
                    "overlap_first_to_second": [list(pair) for pair in item.overlap],
                    "canonical_key_sha256": sha256(key.encode()).hexdigest(),
                })

    all_keys = set(counters["r_plus_after_r_id"]) | set(counters["r_id_after_r_plus"])
    position = {name: index for index, name in enumerate(published_order)}
    rows = []
    for key in all_keys:
        name = _name(key, aliases)
        row = {
            "rule": name,
            "r_plus_after_r_id": counters["r_plus_after_r_id"][key],
            "r_id_after_r_plus": counters["r_id_after_r_plus"][key],
            **convenience_flags(representatives[key], chemistry),
            "canonical_key_sha256": sha256(key.encode()).hexdigest(),
        }
        rows.append(row)
    rows.sort(key=lambda row: (position.get(row["rule"], len(position)), row["rule"]))

    computed_by_name = {row["rule"]: row for row in rows}
    discrepancies, annotation_discrepancies = [], []
    for name in dict.fromkeys([*published_order, *computed_by_name]):
        expected = published.get(name, {})
        actual = computed_by_name.get(name, {})
        for field in ("r_plus_after_r_id", "r_id_after_r_plus"):
            got, want = actual.get(field, 0), expected.get(field, 0)
            if got != want:
                discrepancies.append({"rule": name, "field": field, "published": want, "computed": got})
        if expected and actual:
            for field in ("violates_C1", "violates_C2"):
                if actual[field] != expected[field]:
                    annotation_discrepancies.append({"rule": name, "field": field,
                        "published": expected[field], "computed": actual[field]})

    report = {
        "schema_version": 2,
        "status": "COMPUTED_DIFFERS_FROM_PUBLISHED_TABLE",
        "semantics": {
            "phase": "plain-rule composition before condition completion",
            "application_conditions": "none on the exported row rules",
            "why": (
                "TCS Table 1 classifies plain DPO composites under the declared "
                "chemical filters. Nested/NoEdge conditions enter only in the later "
                "condition-completion calculation."
            ),
        },
        "inputs": {
            "reaction": "aldol_plain.gml",
            "observable": "query_unconditioned.gml",
            "chemistry": "chemistry.json",
        },
        "method": {
            "composition": "MolFoundry common-overlap DPO candidate enumeration",
            "automatic_subset": "edge-reflecting overlaps, computed from each overlap map",
            "chemical_filter": "balanced per-atom valence and declared maximum valence",
            "classification": "condition-sensitive canonical rule isomorphism",
            "reported_constraints": "C1/C2 evaluated from the computed rule interfaces",
        },
        "rows": rows,
        "totals": {
            "computed": {
                order: sum(counter.values()) for order, counter in counters.items()
            },
            "edge_reflecting_automatic": {
                order: sum(counter.values()) for order, counter in automatic.items()
            },
            "published": {"r_plus_after_r_id": 12, "r_id_after_r_plus": 18},
        },
        "non_edge_reflecting_contributions": non_reflecting,
        "published_comparison": {
            "oracle": "article-table.json (literal transcription, retained unchanged)",
            "alias_mapping": _load_json(HERE / "article-manual-aliases.json")["mapping"],
            "normalized_published_rows": list(published.values()),
            "constraint_annotation_discrepancies": annotation_discrepancies,
            "alias_audit_scope": "Plain-rule identities and counts only; C1/C2 annotations are not audited by the alias check.",
            "discrepancies": discrepancies,
            "interpretation": (
                "The historical automatic 9/18 subset is reproduced. Exhaustive current "
                "MolFoundry enumeration gives 13/18 rather than the manually augmented "
                "12/18 publication table. The separate 24 September 2026 audit "
                "reproduces 13/18 with MØD 1.0.0.7 and verifies all 31 contributions "
                "by rule isomorphism and multiplicity; it also gives a direct "
                "graph-edit witness for the additional plain class. See "
                "evidence/aldol-mod-audit-2026-09-24 in the paper checkout. "
                "This runner does not rerun that independent audit."
            ),
        },
    }
    return report, representatives


def compute_article_table() -> dict:
    return _analyze()[0]


def render_computed_table(report: Mapping | None = None) -> str:
    report = dict(report or compute_article_table())
    headings = ("rule", "|r+ * rid|", "|rid * r+|", "C1", "C2")
    body = [
        (
            row["rule"], str(row["r_plus_after_r_id"] or ""),
            str(row["r_id_after_r_plus"] or ""),
            "X" if row["violates_C1"] else "",
            "X" if row["violates_C2"] else "",
        )
        for row in report["rows"]
    ]
    widths = [max(len(headings[i]), *(len(row[i]) for row in body)) for i in range(5)]
    lines = [
        "Computed aldol composition table (published Table 1 comparison)",
        "  ".join(headings[i].ljust(widths[i]) for i in range(5)),
        "  ".join("-" * width for width in widths),
        *("  ".join(row[i].ljust(widths[i]) for i in range(5)) for row in body),
    ]
    totals = report["totals"]
    lines.extend([
        f"Computed totals: {totals['computed']['r_plus_after_r_id']}/"
        f"{totals['computed']['r_id_after_r_plus']}; "
        f"edge-reflecting automatic subset: "
        f"{totals['edge_reflecting_automatic']['r_plus_after_r_id']}/"
        f"{totals['edge_reflecting_automatic']['r_id_after_r_plus']}; "
        "published: 12/18.",
        f"Published comparison: {len(report['published_comparison']['discrepancies'])} "
        "count discrepancies after verified alias merging; C1/C2 annotation differences are separate in the JSON.",
    ])
    return "\n".join(lines)


def export_computed_table(output: Path) -> dict:
    report, representatives = _analyze()
    rows_dir = output / "computed-table" / "rows"
    rows_dir.mkdir(parents=True, exist_ok=True)
    by_digest = {
        sha256(key.encode()).hexdigest(): rule for key, rule in representatives.items()
    }
    for row in report["rows"]:
        relative = Path("rows") / f"{row['rule']}.gml"
        (output / "computed-table" / relative).write_text(
            by_digest[row["canonical_key_sha256"]].getGMLString(), encoding="utf-8"
        )
        row["gml"] = str(relative)
    target = output / "computed-table" / "table.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "computed-table" / "README.md").write_text(
        "# Computed plain composition table\n\n"
        "These GML rules intentionally have **no `constrain` blocks**.  They are "
        "the plain DPO composites classified in TCS Table 1, before application-"
        "condition completion.  `chemistry.json` supplies the external valence and "
        "C1/C2 filters used to accept and annotate rows; it does not mutate a plain "
        "row into a conditional rule.\n\n"
        "For the later condition-sensitive overlap, inspect "
        "`../selected-composites/conditioned/` and "
        "`../selected-commutators.json`.  Those files transport the NoEdge condition "
        "through both composition orders.  This artifact does not claim to export "
        "the paper's complete globally conditioned aldol CTMC.\n", encoding="utf-8"
    )
    return {
        "status": report["status"],
        "rows": len(report["rows"]),
        "table": target.relative_to(output).as_posix(),
        "readme": "computed-table/README.md",
        "gml_directory": rows_dir.relative_to(output).as_posix(),
        "totals": report["totals"],
    }
