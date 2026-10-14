# Growing ring: a bounded search is not a closure theorem

The [literal model](../../examples/structural-nac-oscillator/model.json)
lists the rules, guards, rates, and selected `Q0,Q1,Q2` patterns. Read those
inputs before [`run.py`](run.py). It runs the ordinary composition worklist
fresh and passes that result to
[`structural_nac_common.analyze`](../../src/mfobs/models/structural_nac_common.py).

The `linear` result retains a pending frontier at the stated budget; the
`pjc` result records a selected-coordinate non-lumpability witness but does
**not** assert a finite all-size event quotient or EMGF operator. This is an
honest unknown, not a negative closure proof. `run_with_evidence` is the
same call used by the artifact profile to retain the raw worklist.

```sh
PYTHONPATH=src python -m mfobs.example_cli growing-ring \
  --query-file examples/growing_ring/query.json \
  --out /tmp/growing-ring.json
```

Any connected cyclic-word identity `Q0`, `Q12`, `Q201`, and so on can be
selected with the same `--observable NAME=FILE` interface. The word is read
from the GML and checked against its name; no fixed `Q0,Q1,Q2` table is used.
The expected outcome can still be an exact bounded prefix with pending work.
Omitting `--query-file` selects this same checked-in file as the named
example's default; the explicit spelling above makes the scientific input
visible in the artifact command.
