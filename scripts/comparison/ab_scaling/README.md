# A/B composition scaling benchmark

This benchmark compares MolFoundry 0.1.9 and the FRAGGER graph-rewriting core on the same positive directed A/B path model. The [native measurements](../../../notebooks/visual/evidence/ab-scaling/summary.csv) supply the paper's timing table. They are an implementation comparison on this workload, not an asymptotic-complexity or uniform-speedup claim.

## Workload

For length n, the selected observable is the directed labelled path ABAB... with n vertices. The runner also supports an A^n control family. No host-state enumeration is timed.

Both tools receive six positive terminal-extension rules. For p in {R,A,B} and s in {A,B}, each rule replaces p -> T with p -> s -> T. It preserves p and T, deletes the terminal edge, and adds the A/B vertex and two directed edges. No vertex is deleted, so FRAGGER's SPO and MolFoundry's DPO semantics coincide on this workload.

For each rule, both composition orders are enumerated: twelve calls per observable. FRAGGER uses native `leftUnions` and `rightUnions` on directed multigraphs. MolFoundry uses the public-API typed-incidence adapter with `connected=False`, `maximum=False` and `include_empty=True`. Its incidence constraints enforce only the graph signature, without path-family pruning such as acyclicity or degree bounds.

Raw accepted multiplicities and all twelve per-order counts are recorded before canonical grouping/cancellation. The timed MolFoundry workload includes the adapter's canonical grouping and representative GML serialization. Startup, imports, Scala compilation and JVM startup are excluded and recorded separately. Count agreement alone does not establish equality of all composite expressions.

## Versions and protocol

Both tools were measured serially on native Linux x86_64 on `algochem-compute` (Intel Xeon w7-2595X):

- MolFoundry 0.1.9, Python 3.12.3; measured 14 October 2026.
- FRAGGER `rhz/graph-rewriting` at `0d9595c5d6e8918fc605f7d65778568a74f5e2d4`, Scala 2.13.2, Temurin Java 11.0.32.1 with `-Xmx4g`.

Each size uses a fresh process, four small global warm-ups, two case warm-ups and five timed repetitions. The reported statistic is the median. JSON records retain every timing sample, environment, engine/wrapper/adapter hashes, compiler identities and per-order counts.

The comparison covers every even n=2 through 28. FRAGGER's samples from 9 October 2026 are retained unchanged. All twelve ordered counts agree with the MolFoundry results at each size. The table also includes MolFoundry-only sizes n=48,64,96,128.

## Interpretation

The pinned FRAGGER core constructs intersections by enumerating every subset size of the first graph's vertices and then compatible image choices. The workload exposes this implementation-level scaling cost while keeping the rules and requested compositions fixed. Smaller sizes include cases where FRAGGER is faster.

The result establishes neither a uniform speedup nor a general complexity separation between minimal gluings and DPO overlaps. A timeout of the full process, including warm-ups and all repetitions, is not a per-workload lower bound.

## Reproduction

From the package directory, prepare the pinned FRAGGER source and compiler:

```sh
sh scripts/comparison/ab_scaling/prepare_fragger.sh /tmp/fragger-ab-benchmark
```

To rerun both tools, use a Python environment with MolFoundry 0.1.9 installed and Java 11. Both tools must run on the same native host:

```sh
python scripts/comparison/ab_scaling/run_benchmark.py \
  --fragger-checkout /tmp/fragger-ab-benchmark \
  --java /path/to/java11/bin/java --python /path/to/venv/bin/python \
  --molfoundry-version 0.1.9 --families alt \
  --lengths 2 4 6 8 10 12 14 16 18 20 22 24 26 28 \
  --warmups 2 --repetitions 5 --timeout 14400 --out /tmp/ab-scaling-native
```

The runner defaults to MolFoundry 0.1.9. The installed version must match the requested version. The runner checks the pinned FRAGGER source and compiler and all twelve ordered counts and their totals. It records actual engine, wrapper and adapter hashes.

Optional diagnostic plots:

```sh
python scripts/comparison/ab_scaling/plot_results.py /tmp/ab-scaling-native/results.json
```
