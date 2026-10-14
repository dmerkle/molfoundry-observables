# Model inputs

A model supplies a graph family, rated rules and observable queries. Start
with the [native Python inputs](public-api.md) and
[graph tutorial](../tutorials/graph_closure_tutorial.py).

## File-based manifests

The manifest interface below loads inputs from files for adapters that support
it. `load("examples/eight-slots/model.json").specification()` exposes the same
`omega`, `rules` and `observables` fields as the Python examples. `ModelManifest` stores the file declarations; `ManifestBackend(manifest)`
executes their counting adapter. Numerical settings remain separate from
symbolic synthesis.

Paths in `model.json` are relative to its directory. The
[schema](../schema/model.schema.json) and [loader](../src/mfobs/inputs.py)
check the structure, files, identifiers and rational match factors.

| Field | Meaning |
|---|---|
| `schema_version` | Manifest format, version 1. |
| `id`, `title` | Model identifier and description. |
| `adapter` | Python class `module:Class` supplying the graph family and counting operations. |
| `rules` | Ordered rules with `id`, `file`, `rate`, `per_match` and `kind`. |
| `observables` | Available query GML files with identifiers. |
| `initial_observables` | Selected queries that seed the search. |
| `composition` | Composition backend and overlap options. |
| `forbidden` | Model-specific negative local constraints, or `null`. |
| `state_language` | Graph-family assumptions. |
| `specifications` | Optional chemical signature and domain specifications. |
| `reduction` | Counting procedure and its budgets. |
| `normalisation` | Raw-count normalization for the numerical readout. |
| `initial_mixture`, `initial_state` | Numerical initial conditions. |
| `parameters` | Numerical rate values used by the solver. |
| `reference` | Independent comparison data used by verification. |

The composition mode is `public`, or `restricted` with a `forbidden` file
for graph-family pruning. Both use the installed MolFoundry public API.

A rule's `per_match` factor multiplies its rate for each admissible injective
match. Matches leading to isomorphic products retain their multiplicity.
The `kind` field classifies a rule; it does not prove chemical validity.
Numerical initial conditions do not restrict symbolic synthesis to the size
of that particular mixture. Reference equations are verification inputs.

## Rules, conditions and graph families

For a rule `L <- K -> R`, GML `context` contains the preserved material,
`left` the deleted material and `right` the created material. Thus the full
input is `context + left` and the full output is `context + right`.
An identity rule represents a query pattern.

Rules and queries are partial patterns. A missing edge is forbidden only
when a condition says so. Unshown molecular neighbours must not be filled
in automatically. Adjacency conditions count matching neighbours in the
whole host graph, including neighbours outside the positive pattern.

The allowed graph family is a separate assumption. For example, one ring
or a unique growing path end can justify counting identities. A backend must
establish that its reductions hold throughout that family and that the rules
preserve it. Finite host tests alone do not prove an unbounded-family claim.
The [semantic contract](../src/mfobs/core/semantic_contract.py) records exact,
zero and unsupported reductions together with their justification.

Changing rates changes a numerical experiment. Changing rules, conditions or
the graph family may invalidate its counting argument. See
[other rules and backends](new-model.md) for the available entry points.
