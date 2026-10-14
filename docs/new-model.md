# Other rules and backends

The [graph tutorial](../tutorials/graph_closure_tutorial.py) supplies native
MolFoundry GML rules and a query graph in `ModelInputs`, declares Ω with
`SimpleGraphs`, and selects `SimpleGraphBackend`. It calls Algorithms
1 and 2 separately:

```sh
PYTHONPATH=src:. python -m tutorials.graph_closure_tutorial
```

For a packaged family, `backend.specification()` exposes its Ω, rated MolFoundry
rules and observable expressions; see the [Python input API](public-api.md).
Its counting backend accepts only supported inputs and rejects other edits.
For example, changing `MOTIF` in the
[ABA notebook](../notebooks/examples/02_aba.ipynb) selects another query
while keeping the rules and counting backend fixed.

## Change a bundled example

Start with [ABA](../tutorials/aba_analysis.py): change `MOTIF` to `"ABBA"`
while keeping Ω and the path rules fixed. [Molecular ROP](../tutorials/rop_analysis.py) exposes `MOTIF` and `GATED`.
Changing `"ABA"` to `"ABBA"` reuses the same molecular rules and counting
backend. Its [inputs](../src/mfobs/models/molecular_rop/inputs.py) contain the
monomer, chain and reaction constructors; its `counting.py` and `events.py`
justify the reductions on that constructor family.

[Chlorination](../tutorials/chlorination_analysis.py) compares the same rules
with and without the neighbour NAC. Its growth factor is `1/(3!·4!)`, cancelling
labelled-hydrogen match permutations. Its backend supports those two rules
and the declared four-site family; arbitrary edited chemistry is not covered.

To introduce a different graph family, supply Ω, rules, queries and the needed
exact counting backend. The generic worklists and arithmetic remain unchanged.
Native MolFoundry matching and composition alone do not supply identities on Ω.

## Inspect one commutator

For literal rule and query GML, the CLI composes both orders and exports the
signed terms:

```sh
python run.py derive \
  --observable examples/gated-path/observables/vertices.gml \
  --rule examples/gated-path/rules/subdivide.gml \
  --out out/one-commutator
```

This operation needs no family-specific count reducer. Its output is a
commutator calculation; additional counting arguments are needed to turn
those terms into observable equations.

## Supply a counting backend

The [algorithm map](../notebooks/visual/CORE.md) identifies the common linear
and support searches and their concrete backend interfaces. A backend
accepts the declared graph family and supplies generator images, exact
arithmetic and counting identities. `CountingBackend.activity(rule)` and
`CountingBackend.image(pattern, rule)` feed the shared support worklist. It must report unsupported reductions explicitly. Event
extraction additionally needs a justification for the complete distribution
of coordinate jumps. The composition-based Algorithm 1 adapter also uses
`counter.reduce(rule)` and `counter.compose(pattern, rule)`; the
[ABA counting module](../src/mfobs/models/aba/counting.py) shows these operations.

For the optional [file-based manifest interface](inputs.md),
[`DynamicBasisAdapter`](../src/mfobs/core/synthesis.py) connects graph
composition and certified reduction to the shared linear loop.
[`SemanticContract`](../src/mfobs/core/semantic_contract.py) records the
reduction scope and rule-preservation arguments. The
[Eight sites adapter](../src/mfobs/models/eight_slots.py) is a concrete example.
Its invariants belong to that model and cannot be assumed for arbitrary GML.

The canonical identity adapter performs reduction only by graph-rule identity:

```sh
python run.py run path/to/model.json --max-expansions 20 --out out/model
```

Use `"adapter": "mfobs.core.identity:Adapter"` in that manifest. An unfinished
search reports `UNKNOWN` with pending work. For numerical integration, supply
initial expectations for every resulting coordinate. A finite exact system
and an all-state counting argument are separate requirements from passing
finite implementation tests.
