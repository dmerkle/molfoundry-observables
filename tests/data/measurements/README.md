# Native measurements

`query-discovery.json` records three fresh serial processes per motif on the prepared software source, including source/runtime hashes and independent checks. `composition.json` retains the measured MolFoundry/FRAGGER comparison; the timed composition source files are unchanged. Full analysis outputs are reproducible with the [query](../../../scripts/comparison/query_reuse/README.md) and [composition](../../../scripts/comparison/ab_scaling/README.md) runners. These measurements are evidence, not synthesis inputs.

The fifteen raw query results are compressed under `query-runs/`; their checksums are recorded in the summary. `fragger-original.json` is the unchanged source record identified by the composition summary.

`source-snapshot.tar.gz` retains the exact source files and measurement runner
whose hashes are recorded in `query-discovery.json`. This evidence remains
verifiable independently of the public Git history. The maintained sources
relocate the text formatter and use POSIX manifest paths; measured algorithms
and counting code remain unchanged. `scripts/comparison/query_reuse/records.py` checks the
snapshot hashes and the maintained analysis sources, allowing only those
bookkeeping changes. The timing samples and raw analyses
have not been rewritten or remeasured.
