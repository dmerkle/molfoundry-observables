"""Run one example analysis; each example keeps its own readable input and call.

From the implementation root: PYTHONPATH=src python -m mfobs.example_cli EXAMPLE --out FILE
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from mfobs.query_input import from_arguments, from_file
from examples.brusselator.run import run as brusselator
from examples.chlorination.run import run as chlorination
from examples.growing_ring.run import run as growing_ring
from examples.remy.run import run as remy
from examples.ring8.run import run as ring8
from examples.aba.run import run as aba
from examples.molecular_rop.run import run as molecular_rop
from examples.tcs.replay import run as tcs
from mfobs.formatting import render


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_QUERY_ARGUMENTS = {
    "brusselator": (
        "X=examples/brusselator-growth/observables/X.gml",
        "Y=examples/brusselator-growth/observables/Y.gml",
    ),
    "ring8": ("P0=examples/ring8/observables/occupied.gml",),
}
DEFAULT_QUERY_FILES = {
    "growing-ring": ROOT / "examples/growing_ring/query.json",
    "aba": ROOT / "examples/aba/query-gml.json",
    "remy": ROOT / "examples/remy/query.json",
}


def _default_query_record(example: str):
    if example in DEFAULT_QUERY_FILES:
        return from_file(DEFAULT_QUERY_FILES[example])
    values = DEFAULT_QUERY_ARGUMENTS.get(example)
    return from_arguments(values, base=ROOT) if values else None


def _tcs(*, queries=None):
    if queries is not None:
        raise ValueError("the local TCS diagnostic has no selected ODE observable")
    return tcs()


EXAMPLE_RUNS = {
    "tcs": _tcs,
    "brusselator": brusselator,
    "chlorination-guarded": lambda *, queries=None: chlorination("guarded", queries=queries),
    "chlorination-ablated": lambda *, queries=None: chlorination("ablated", queries=queries),
    "ring8": ring8,
    "aba": aba,
    "rop": molecular_rop,
    "remy": remy,
    "growing-ring": growing_ring,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("example", choices=EXAMPLE_RUNS)
    parser.add_argument("--out", required=True, type=Path)
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument(
        "--observable", action="append", default=[], metavar="NAME=FILE[@COEFFICIENT]",
        help=(
            "replace the example's selected query; repeat the same NAME for a "
            "finite expression such as Rémy's three stem components"
        ),
    )
    query_group.add_argument(
        "--query-file", type=Path,
        help="JSON form of the same observable selection",
    )
    args = parser.parse_args()
    if args.out.exists():
        parser.error("output exists; choose a new path")
    try:
        explicit_queries = (
            from_file(args.query_file) if args.query_file else
            from_arguments(args.observable) if args.observable else None
        )
        queries = explicit_queries or _default_query_record(args.example)
        result = EXAMPLE_RUNS[args.example](queries=queries)
    except ValueError as exc:
        parser.error(str(exc))
    if queries is not None:
        result["query_input"] = queries.record()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(render(result, name=args.example))
    print(f"\nFull derivation record: {args.out}")


if __name__ == "__main__":
    main()
