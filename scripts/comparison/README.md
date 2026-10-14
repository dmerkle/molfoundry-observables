# FRAGGER comparisons

The paper's external comparison is the two-legged DNA walker from Danos et al.
Both tools derive the extended/compressed count equations. The retained
[evidence](../../notebooks/visual/evidence/dna-walker/README.md) records the inputs,
normalization and comparison with the published velocity.

## ROP population projections

`RopPopulationHarness.scala` compares two edgeless population projections:

| Preset | Rules | Query | FRAGGER basis |
|---|---|---|---|
| `baseline` | `empty -> M`, `M -> N` | `N` | `N, M, 1` |
| `noedge` | `empty -> M`, `O -> C`, `C -> O`, `O+M -> O+N` | `N` | `N, M+O, O, M, 1` |

The correspondence is `N -> N_B`, `M -> M_B`, `O -> X`, `C -> 1-X`.
The gated population language supplies `O+C=1`, `O²=O`, `C²=C`, `OC=0`.
The runner compares each FRAGGER image with MolFoundry's edgeless compositions
and the molecular ROP event generator, preserving raw injective multiplicities.
These are population projections; they do not assert identical molecular input
or give FRAGGER native molecular NAC semantics.

Use Java 11 and the pinned upstream checkout:

```sh
git clone https://github.com/rhz/graph-rewriting.git /tmp/fragger-rop
git -C /tmp/fragger-rop checkout --detach 0d9595c5d6e8918fc605f7d65778568a74f5e2d4
mkdir -p /tmp/fragger-rop/compiler
for part in compiler library reflect; do
  curl -fL "https://repo.maven.apache.org/maven2/org/scala-lang/scala-$part/2.13.2/scala-$part-2.13.2.jar" \
    -o "/tmp/fragger-rop/compiler/scala-$part-2.13.2.jar"
done
PYTHONPATH=src python scripts/comparison/run_rop_population_comparison.py \
  --checkout /tmp/fragger-rop --java /path/to/java11/bin/java --out /tmp/rop-population
```

The runner checks the upstream commit, clean core sources and compiler JAR
hashes. `comparison.json` contains inputs, mapped images, subprocess outcomes
and provenance. Its phase timings are diagnostics, not the paper's measurements.

## Bounded directed-word comparison

`RopSequenceHarness.scala` uses one directed `R -> word -> T` chain and isolated
A/B pool vertices. It supplies feed/append rules and seven count queries.
The constant edge label is omitted; direction, endpoints, vertex labels and rule
interfaces are preserved. The pinned core mishandles deletion of labelled edges,
so the labelled version is retained only as a diagnostic.

The finite input check covers 135 word/pool states. Bounded FRAGGER runs produce
12 and 24 equations; 7 and 8 respectively map to supported word/suffix counts
and agree exactly with MolFoundry. The other equations contain unsupported
branching or multiple word components and remain `UNMAPPED`. The overall result
is `COMPARISON_INCOMPLETE`; a capped search proves neither closure nor nonclosure.

```sh
PYTHONPATH=src python scripts/comparison/run_rop_sequence_comparison.py \
  --checkout /tmp/fragger-rop --java /path/to/java11/bin/java --out /tmp/rop-sequence
PYTHONPATH=src python -m pytest -q tests/test_rop_population_comparison.py \
  tests/test_rop_sequence_comparison.py
```
