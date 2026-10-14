# Reproduction profiles

`paper_evidence_manifest.json` declares the model settings, process budgets and expected outcomes. `reproduce_paper.py` executes each analysis in a fresh process and writes its result, checks and logs into the specified directory.

```sh
python scripts/reproduction/reproduce_paper.py run --profile smoke --output out/smoke
python scripts/reproduction/reproduce_paper.py run --profile paper --output out/paper
```

Use a new or empty output directory. `run.json` records every stage and its outcome; `claims.md` links the individual results. Expected UNKNOWN statuses count as reproduced only when their structured checks pass. They do not establish non-closure.

The [paper claim map](../../docs/paper-claims.md) describes the scientific scope. The root [Dockerfile](../../Dockerfile) packages the pinned runtime and the same profile. `supplementary` adds fixed Ring-8, complete-suffix and external comparison workloads; those require the pinned FRAGGER checkout and Java options shown by `--help`.

Recorded timings belong to their recorded source revision and machine. Local process durations are reproduction diagnostics.
