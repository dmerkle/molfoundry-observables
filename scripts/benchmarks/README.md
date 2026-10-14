# ROP query scaling

`rop_scaling.py` measures a fixed complete-suffix refinement of the ROP model. It keeps the literal A/B chemistry and feeds fixed while adding
nested motif queries at terminal context depths zero through three, each without
and with the existing gate. Every observation is a fresh serial Python process,
compiles the same refined rules and calls the shared linear and jump analysis
with linear budget 24. The parent enforces a process timeout. Process failures remain rows.

Smoke profile (letters, both gate settings):

```sh
PYTHONPATH=src python scripts/benchmarks/rop_scaling.py run \
  --profile smoke --repetitions 1 --timeout 600 --output /tmp/rop-scaling-smoke
```

Full eight-case matrix:

```sh
PYTHONPATH=src python scripts/benchmarks/rop_scaling.py run \
  --profile full --repetitions 1 --timeout 600 --output /tmp/rop-scaling-full
```

The `abba` query includes `ABA` and `ABBA`; it is not the older
dyads-plus-ABBA timing diagnostic. A later frozen artifact may use
`--repetitions 3` on the identical matrix. Single-run values are diagnostic
measurements without error bars, speedup, or statistical claims.

## Match counting

`python scripts/benchmark_molfoundry_counts.py --carbons 24 54 104` compares
MolFoundry's direct count with DG application for a one-carbon query. Each
measurement runs in a fresh process and checks the returned count.
