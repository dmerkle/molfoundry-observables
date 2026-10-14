# Local slot experiment

These inputs reuse the original Ring8 rule GMLs and occupied-slot seed, but
**not** the finite Ring8 scaffold, state enumeration, orbit vectors or reducer.

- `model-local.json`: only the declared ownership/type/capacity axioms.
- `model-degree2.json`: additionally a maximum next-degree of two; still no
  fixed cycle length, connectedness, or requirement to be a cycle.
- `model-no-coop.json`: basal binding and unbinding only, a closed positive control.

The invariant JSONs are trusted supplied assumptions. The production reducer
implements this particular fragment soundly but incompletely; it does not
discover the assumptions. Pattern-size and expansion budgets mean UNKNOWN,
not a semantic bound on hosts or a truncated closed ODE.

See [the experiment and proof](../../docs/local-slot-experiment.md).

Use the environment from the [package setup instructions](../../README.md#install-and-open) and run these commands from the package directory.

```sh
python scripts/run_local_slots.py --model all --max-expansions 7
python scripts/check_local_slots.py
```
