#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    data = json.loads(args.results.read_text(encoding="utf-8"))
    out = (args.out or args.results.parent).resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows = data["comparison"]
    for family in ("alt", "same"):
        subset = [row for row in rows if row["family"] == family]
        if not subset:
            continue
        subset.sort(key=lambda row: row["n"])
        x = [row["n"] for row in subset]

        fig, ax = plt.subplots()
        ax.plot(
            x,
            [row["molfoundry_median_seconds"] for row in subset],
            marker="o",
            label="MolFoundry",
        )
        ax.plot(
            x,
            [row["fragger_median_seconds"] for row in subset],
            marker="s",
            label="FRAGGER",
        )
        ax.set_yscale("log")
        ax.set_xlabel("Observable length n")
        ax.set_ylabel("Median composition time (s)")
        ax.legend()
        ax.grid(True, which="both", axis="y", linewidth=0.5)
        fig.tight_layout()
        stem = out / f"ab-composition-runtime-{family}"
        fig.savefig(stem.with_suffix(".pdf"))
        fig.savefig(stem.with_suffix(".png"), dpi=300)
        plt.close(fig)

if __name__ == "__main__":
    main()
