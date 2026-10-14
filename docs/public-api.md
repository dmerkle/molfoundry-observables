# Public API and backend boundary

There are three distinct parts:

- **Model inputs:** Ω, native MolFoundry rules with rates, and observable patterns.
- **Counting backend:** exact identities on that Ω and, for events, its jump argument.
- **Shared framework:** input types in `mfobs.api` and Algorithms 1/2 in `mfobs.core`.

`mfobs.models.<example>.inputs` defines a bundled model. Its `counting` and
`events` modules supply family-specific mathematics; they are not the generic
algorithms. `mfobs.applications` contains presets that assemble those parts.
The [algorithm map](../notebooks/visual/CORE.md) links the actual files.

Start with the [graph tutorial](../tutorials/graph_closure_tutorial.py), which
writes Ω, rule GML and query GML directly before calling either algorithm.

`ModelInputs` contains only Ω, rules and observables; it has no execution method.
`Omega` declares the admissible state family. Rules are native MolFoundry rules
with rate channels and per-match factors. Observables are named graph counts
or rational sums of identity-rule counts; `identity_terms()` exposes their
coefficients and native rules. Rémy's three patterns form one selected query.
Brusselator molecule counts include their symmetry normalization.
`rule_gmls()` and `query_gmls()` provide the pattern drawings.

Each model declares Ω in `state_family()` beside its constructors; the short
[`application specifications`](../src/mfobs/applications.py) pair them with rules
and selected queries. Eight sites and the ring examples read their existing
`state_language` manifest input. [`state_space.py`](../src/mfobs/state_space.py)
provides the common type, manifest reader and chemical restrictions. Chemical families share
[`ordinary-covalent-v1.json`](../src/mfobs/data/ordinary-covalent-v1.json), with
explicit coordination and boundary-marker extensions where needed.

The same interface covers `ABA`, `MolecularROP`, `RemyTrees`, `EightSlots`, `Chlorination`, `Ring8`,
`BrusselatorROP`, `GrowingRing` and `DNAWalker`. Unsupported edits to Ω, rules or
queries are rejected before analysis. A family declaration alone does not
supply a counting theorem: each backend implements its stated reductions.
Discovered support belongs to the result, not these inputs.

## Define a model

The [graph tutorial](../tutorials/graph_closure_tutorial.py) is a complete
runnable example with inline MolFoundry GML for the rule and query graph.
`RatedRule` associates the parsed MolFoundry rule with a rate; `Observable`
names the count of a MolFoundry graph pattern. `SimpleGraphs` declares the
allowed graph family. `SimpleGraphBackend(Omega)` supplies its counting operations.
Its `linear_closure(inputs)` and `jump_closure(inputs)` methods run Algorithms 1
and 2 separately; `analyze(inputs, backend=backend)` runs both.

The public `analyze` call returns `AnalysisResult` with separate `linear`, `support` and `pjc`
outcomes. A missing outcome is `None`; support closure does not imply event or
linear closure. Each computed outcome exposes `status` and `data`. `record()`
exports the full numerical/evidence record consumed by the drawings and solvers.

File-based inputs load as `ModelManifest`. Its `specification()` returns the
same native `ModelInputs`; `ManifestBackend(manifest)` runs the declared counting
adapter. Loading a manifest does not run an analysis.

A new graph family may need its own exact count reductions or event partitions.
A state-membership test alone cannot supply those identities. See the
[new-model guide](new-model.md).
The [algorithm map](../notebooks/visual/CORE.md) links the shared worklists and
the operations supplied by each backend.

## MolFoundry boundary

The package uses the public Python API of MolFoundry 0.1.9:

| Operation | MolFoundry API |
|---|---|
| Spans and conditions | Rule vertices/edges, constraints and `matchCondition` |
| Canonical grouping | `canonKey` |
| Identity observables | `Rule.identity` |
| Composition | `RCEvaluator.iterate(..., onlyUnique=False)` |
| Overlap maps | `CompositionVertex.firstId` / `secondId` |
| Jump projection | `Rule.dpoJumpClosure`, `Rule.sqpoJumpClosure` |
| Count or enumerate matches | `Rule.numMatches`, `Rule.matches` |
| Apply a supplied match | `DG.Builder.applyMatch` |
| Rule drawings | `Rule.svg` / `Rule.print` with `GraphPrinter` |

Canonical grouping retains raw multiplicity. Local graph-family restrictions
are compiled into `IncidenceConstraint`s so MolFoundry can prune overlaps.
Nonlocal restrictions, such as directed cycles, use public filter predicates.
The `maximumCandidateOverlaps` limit counts intermediate graphs passed to the
predicate; it does not bound work inside one iterator step.

Concrete-state counting uses `ApplyMode.UseAny`: only the supplied graph
copies are available. Conditions, injectivity and symmetry multiplicities are
handled by MolFoundry. An empty pattern has one match. `DG.HyperEdge.numMatches`
reports match multiplicities for a DG edge; `vertexMaps` gives atom correspondences for products and
highlighting, not the rule matches themselves.

Graphs and rules are MolFoundry objects throughout the analysis.
[`core/graph.py`](../src/mfobs/core/graph.py) provides construction and labelled-ID
access helpers; [`core/rules.py`](../src/mfobs/core/rules.py) reads native spans.
[`core/molfoundry_execution.py`](../src/mfobs/core/molfoundry_execution.py)
delegates matching and application to MolFoundry and preserves host vertex IDs
in products. NetworkX is used for structural graph algorithms. MolFoundry
provides SMILES and ordinary molecule views; RDKit supplies exact bond highlights,
selected explicit hydrogens and charged aromatic views, and independent chemical checks.
