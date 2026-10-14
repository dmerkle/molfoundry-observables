# Validation

Synthesis reads rules, queries and the state-language contract. Reference
results are used only by the separate verification step.

- Composition checks compare canonical rule multiplicities and condition
  transport against direct applications on concrete graphs.
- Fixed-ring checks enumerate all 256 states and compare rate-separated
  generators with an independent bit-flip model.
- Chain-family checks compare structural count formulas with finite host
  matching. Finite checks support, but do not replace, the all-size
  [counting and preservation arguments](reusable-family.md).
- Chlorination NAC checks compare joint jumps and rates with direct MolFoundry
  applications and an independent concrete reference.

```sh
bash scripts/reproduce.sh
```

Unsupported reductions or exhausted analysis budgets return `UNKNOWN`, not a
closure certificate. Input and verification failures are errors. See the
[public API boundary](public-api.md) for remaining engine limitations.
