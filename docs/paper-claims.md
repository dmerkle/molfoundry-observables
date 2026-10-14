# Claims and checks

The [CLI guide](../tutorials/README.md) gives the shortest path from inputs to
results; the [algorithm map](../notebooks/visual/CORE.md) identifies Algorithms 1
and 2 and their family-specific counting/event arguments.

| Example | Inspect | Independent checks |
|---|---|---|
| Conditional composition | [TCS inputs](../examples/tcs/README.md) | `tests/test_aldol_current_mod.py`, `examples/tcs/replay.py`; literal host applications and signed composition contributions. |
| ABA and query variation | [ABA inputs](../src/mfobs/models/aba/inputs.py) | `mfobs.checks.rop_projection`; bounded word reference and direct graph applications. Algorithm 1 reaches its budget; Algorithm 2 closes. |
| Rémy trees | [Counting argument](remy-trees.md) | `tests/test_remy_trees.py`; recursive tree counts and direct applications. |
| Eight sites | [State family and reduction](eight-slot-topologies.md) | `scripts/check_eight_slots.py`; all 768 occupancy assignments, after symbolic synthesis. |
| DNA walker | [Original rules and recovered exports](../examples/dna-walker/README.md) | `mfobs.checks.dna_walker`; fresh equations and retained FRAGGER calculations. |
| Molecular ROP | [Counting argument](molecular-rop-discovery.md) | `mfobs.checks.molecular_rop`; literal molecular matches, successors and joint jump/rate distributions. |
| Chlorination | [NAC and ablation](../examples/chlorination/README.md) | `tests/test_context_nac_chlorination.py`, `tests/test_context_nac_direct.py`; the two variants differ only in the neighbour NAC. |
| Composition measurements | [Protocol and external setup](../scripts/comparison/ab_scaling/README.md) | Retained samples, source identity and twelve ordered counts at each size. |

Run the semantic suite with the `paper` reproduction profile documented in the
[README](../README.md#check-and-reproduce). It computes analyses from inputs;
large timing repetitions and FRAGGER execution are separate external workloads.
The [native measurements](../tests/data/measurements/README.md) retain source identity, versions, samples and hardware. The [query runner](../scripts/comparison/query_reuse/README.md) reproduces the five-query timing table.

Exactness is relative to the declared graph family, exact reduction identities
and justified event construction. Bounded tests check implementations; they do
not replace the all-family counting arguments. UNKNOWN denotes unfinished or
unsupported analysis, not proof of non-closure. Approximation curves in the two
further examples are identified separately from exact generated results.
