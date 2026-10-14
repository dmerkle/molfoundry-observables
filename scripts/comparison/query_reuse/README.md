# Query discovery measurements

Run five ABA motifs with unchanged rules and counting backend. Each repetition uses a fresh serial process. The reported Algorithm 2 time is support discovery plus event extraction; Algorithm 1 and independent word checks follow the timed call.

From a clean checkout, on the same native host as the composition comparison:

```sh
python scripts/comparison/query_reuse/measure.py --repetitions 3 \
  --cpu-description 'your CPU model' --out /tmp/query-discovery
```

`summary.json` contains all samples, medians, source/runtime hashes and check outcomes. Individual analyses remain in that output directory. Docker/emulation durations are not the paper's native measurements.
