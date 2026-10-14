"""The paper's aldol question: compose the rule with an identity query both ways.

Only the common-overlap terms are needed for the selected local r180 claim.
The paper's broader calculation also considers parallel terms and global Chem
constraints; this example does not claim to reproduce that full CTMC.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import molfoundry as mf

from . import _audit


HERE = Path(__file__).resolve().parent


def load_rules(*, omega=None) -> dict[str, str]:
    """The reaction and identity queries, checked against shared local chemistry."""
    from mfobs.state_space import covalent_compositions

    omega = (omega or covalent_compositions()).require(covalent_compositions())
    inputs = {
        name: (HERE / f"{name}.gml").read_text(encoding="utf-8")
        for name in ("aldol", "query_unconditioned", "query_conditioned")
    }
    for gml in inputs.values():
        omega.validate_rule(mf.Rule.fromGMLString(gml, add=False))
    return inputs


def compose(query, rule) -> dict:
    """Return both public composition orders, retaining overlap multiplicity."""
    common = mf.rcCommon(connected=False, includeEmpty=False)
    evaluator = mf.RCEvaluator([])
    return {
        "before": evaluator.eval(query * common * rule, onlyUnique=False),
        "after": evaluator.eval(rule * common * query, onlyUnique=False),
    }


def run_with_products(*, omega=None) -> tuple[dict, dict]:
    """Compute the checked record and retain products for optional GML export."""
    inputs = load_rules(omega=omega)
    aldol = mf.Rule.fromGMLString(inputs["aldol"], add=False)
    queries = {
        name: mf.Rule.fromGMLString(inputs[f"query_{name}"], add=False)
        for name in ("unconditioned", "conditioned")
    }
    products = {name: compose(query, aldol) for name, query in queries.items()}
    return _audit.verify(inputs, aldol, queries, products), products


def run(*, omega=None) -> dict:
    result, _products = run_with_products(omega=omega)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--force", action="store_true", help="Replace the output record")
    parser.add_argument(
        "--article-table", action="store_true",
        help="Recompute the aldol composition table and compare it with TCS Table 1",
    )
    parser.add_argument(
        "--gml-dir", type=Path,
        help="Also export computed table rows and replay evidence as GML files",
    )
    parser.add_argument(
        "--gml-scope", choices=("selected", "all"), default="selected",
        help="With --gml-dir, optionally export every locally computed composition class",
    )
    args = parser.parse_args()
    if args.out.exists() and not args.force:
        parser.error("output exists; choose a new path or pass --force")
    result, products = run_with_products()
    result["reproducer_sha256"] = sha256(
        Path(__file__).read_bytes() + Path(_audit.__file__).read_bytes()
    ).hexdigest()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PASS: conditioned overlap coefficients −1,+1; separating host counts 0,1")
    print(args.out)
    if args.article_table or args.gml_dir:
        from .inspect import export_gml, render_computed_table
        if args.article_table:
            print("\n" + render_computed_table())
        if args.gml_dir:
            export_gml(
                result, args.gml_dir,
                products=products if args.gml_scope == "all" else None,
                force=args.force,
            )
            print(f"GML export ({args.gml_scope}): {args.gml_dir}")


if __name__ == "__main__":
    main()
