# ABA on growing paths

Run the [ABA analysis](../tutorials/aba_analysis.py) from the package directory:

```sh
python -m tutorials aba
python -m tutorials aba --out /tmp/aba --details
```

Change `MOTIF` or `GATED` in the Python example, or use the
[notebook](../notebooks/examples/02_aba.ipynb).
The rules feed free A/B tokens and append at the unique path end.

From the ABA query, Algorithm 2 discovers `C_ABA, M_A, M_B, S_AB, S_A`,
with `X` added by the gate. Activities introduce the pools; the query's
image introduces `S_AB`, whose image introduces `S_A`. These are discovered
support counts. Event extraction gives eight events without the gate and
ten with it. Algorithm 1 runs independently by composition and count
reduction; both linear searches remain unresolved at dimension 24.

## Counting and events

The unique terminal and interchangeable free tokens give a common coordinate
jump for every match of a source rule at a fixed state. The backend checks
the source rules before using that argument.

Suffix indicators are idempotent. Compatible suffixes multiply to the
longer suffix; incompatible ones multiply to zero, such as `S_A*S_AB=0`.
These are identities of graph counts, not independence assumptions.
The resulting EMGF operator has order two, or three with the gate.

`--out` saves the readable summary, definitions and event table. `--details`
also saves the full record. Support, linear closure and pJC have separate
statuses; exhausted budgets or unsupported reductions retain pending work
and do not certify a complete event model.

[Counting backend](../src/mfobs/models/aba/counting.py) ·
[Event argument](../src/mfobs/models/aba/events.py) ·
[Algorithm map](../notebooks/visual/CORE.md) ·
[Atom-resolved ROP](molecular-rop-discovery.md).
