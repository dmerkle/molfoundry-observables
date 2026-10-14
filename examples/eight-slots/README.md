# Eight sites: symbolic synthesis

The input `examples/eight-slots/model.json` specifies three conditional rules,
the occupied-site query `P0`, and ownership, total-size-eight and degree-two
invariants. It supplies neither a host graph nor a state table. The possible
core topologies are 8, 4+4 and 3+5.

## Run and inspect

From the repository root, with the analysis environment installed:

```sh
python -m tutorials eight-sites --out /tmp/eight-sites
```

Start with `summary.txt`, then `observables.txt` and `ode.txt`. Add `--details`
to retain complete derivation records and literal patterns. The independent
checks run after synthesis; they do not supply its observable basis.

For the complete synthesis and exhaustive occupancy check:

```sh
python scripts/run_eight_slots.py --out out/eight-slots
python scripts/check_eight_slots.py --out out/eight-slots
```

The counting and ownership arguments are documented in
[the derivation guide](../../docs/eight-slot-topologies.md). Optional old full
result archives are retained in the original development repository; the
public software computes its analysis from the declared inputs.
