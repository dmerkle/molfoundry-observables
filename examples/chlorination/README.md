# Chlorination with a neighbour NAC

[`Chlorination()`](../../src/mfobs/applications.py) selects the single-site rule with a genuine `not exists` condition: no neighbouring ring carbon may carry chlorine. `Chlorination("ablated")` removes this condition. Both use the same unbounded tail family with total chlorine four and select only `Z`; support counts and event classes are discovered.

Inspect `rule_gmls()` and `query_gmls()` before calling `analyze(model.specification(), backend=model)`. The [notebook](../../notebooks/examples/03_chlorination.ipynb) compares both rules and their generated dynamics.

```sh
PYTHONPATH=src python -m mfobs.example_cli chlorination-guarded --out /tmp/chlor-guarded.json
PYTHONPATH=src python -m mfobs.example_cli chlorination-ablated --out /tmp/chlor-ablated.json
```

The backend currently selects `Z` only; custom `--observable` selections are unsupported.
