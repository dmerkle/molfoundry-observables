# Aldol addition: conditional composition

This example revisits Behr et al., *Rewriting Theory for the Life Sciences*,
§8.3, DOI [10.1016/j.tcs.2021.07.026](https://doi.org/10.1016/j.tcs.2021.07.026).
The observable counts a molecular pattern whose two designated carbon atoms
must not be bonded. The paper models this condition mathematically and
analyses its effect manually. Here an identity rule with `NoEdge(0,5)`
represents the observable, and MolFoundry computes the conditional compositions.

## Inputs and matched comparison

- [`aldol.gml`](aldol.gml): reaction rule with its own `NoEdge(0,5)`.
- [`query_unconditioned.gml`](query_unconditioned.gml): identity rule without an explicit NoEdge constraint.
- [`query_conditioned.gml`](query_conditioned.gml): the same identity rule with `NoEdge(0,5)`.
- [`chemistry.json`](chemistry.json): balanced per-atom valence and maximum valences C=4, H=1, O=2. The C1/C2 convenience predicates are reported in the published-table comparison, not applied as additional filters here.

[`compare_conditions.py`](compare_conditions.py) composes all common overlaps,
including the empty overlap, transports conditions and applies the same
valence filter in every run. The reaction's constraint stays fixed in the
two main variants; only the identity rule's NoEdge constraint changes.

| Identity rule | Reaction then identity (+) | Identity then reaction (−) | Nonzero signed composition terms |
|---|---:|---:|---:|
| Without explicit NoEdge | 13 | 18 | 13 |
| With NoEdge | 10 | 18 | 18 |

The constraint excludes three reaction-then-identity contributions (plain
classes `r68`, `r73`, `r120`). It also prevents some cancellations: isomorphic
spans with different transported constraints represent distinct conditional
rules. The 18 surviving terms are **composition terms before observable
projection**, not 18 observables or ODEs. A third control, using
[`aldol_plain.gml`](aldol_plain.gml) and the unconditioned identity rule,
reproduces 13/18 and 13 surviving terms.

From the analysis package directory:

```sh
PYTHONPATH=src python -m examples.tcs.compare_conditions --out-dir /tmp/aldol-conditions
```

The output contains `report.json`, literal input copies with hashes, overlap
maps for each accepted contribution, and a GML representative for each
conditional class. Published rule identifiers label the results only; they do not
construct compositions.

## Selected pair and independent host check

[`replay.py`](replay.py) and [`_audit.py`](_audit.py) check the selected
published `r180` pair. Without the identity rule's constraint, the terms
cancel. With it, they have different conditional canonical keys and
coefficients −1,+1. Independently implemented NetworkX matching against the
16-atom host constructed by [`witness_host`](_audit.py#L86) gives counts 0,1 with
the constraint and 1,1 without it. The host is a valence-complete graph, not
a rule span. The overlap maps and host are explicit witnesses for this pair.

```sh
PYTHONPATH=src python -m examples.tcs.replay \
  --out /tmp/tcs.json --article-table --gml-dir /tmp/tcs-gml
```

Inspect `variants.conditioned.selected_composites`,
`selected_net_commutator_coefficients` and `host_witness.input_match_counts`.
The replay command exports the host as `host-witness.gml` in its `--gml-dir`,
together with the selected rule GMLs.
The replay's raw all-composition counts use nonempty overlaps without the
valence filter; they must not be compared directly with the filtered 13/18
and 10/18 totals above.

## Independent MØD comparison

The published table is a comparison oracle, not an input to composition.
The [visual notebook](../../notebooks/examples/01_aldol.ipynb)
computes MolFoundry results afresh and compares individual rule classes and
per-order multiplicities with recorded MØD 1.0.1.100 exports. Both plain and
conditional tables link to GMLs. MØD itself is not rerun by the notebook.

The [25 September replay](reference/mod-1.0.1.100-2026-09-25/README.md) gives **13/18**
plain contributions with `UseBoostCommonSubgraph=False`; all 31 contributions
agree by rule class and per-order multiplicity. The same MØD version gives
9/18 with the flag set to `True`. Thus the overlap implementation, not the
version alone, selects the restricted versus complete enumeration.

The TCS article restricts the computed overlaps to vertex-induced
overlaps and adds three overlaps manually, giving 12/18. Those manual rules
are rA=r73, rB=r120 and rC=r68 (plain-rule isomorphisms; Eq. 104).
After merging these aliases, one additional class remains relative to the
published table: MØD's `r978`, isomorphic to MolFoundry's
`new-f7c8786fa6b1`, also verified using MØD's `Rule.isomorphism`.
Both implementations agree on this control.

The second MØD calculation encodes absence using ordinary `:` edges of valence
weight zero. Translating these exported edges into explicit input-side NoEdge
constraints gives **23 conditional classes**, **10/18 contributions** and
**18 nonzero signed terms**, all matching native MolFoundry composition.
The notebook's conditional tables are self-contained MolFoundry tables: each
numbered row links only its exact native MolFoundry GML, alongside After, Before
and After − before counts. All classes, including cancellations, and totals
remain visible. The original encoded and translated GMLs remain in the separate
archived cross-tool audit. Its correspondence is checked for this example's
resulting rules, not claimed as a general equivalence theorem for absence
encodings or native NAC support in MØD. The notebook retains this check using
[`current_mod.py`](current_mod.py) and the archived GMLs.

`computed-table/table.json` and `computed-table/rows/` from the replay command
provide the plain table and one parseable GML per class. The selected
conditional GMLs and their signed coefficients are exported separately.
To inspect an existing replay without recomputing composition:

```sh
PYTHONPATH=src python -m examples.tcs.inspect /tmp/tcs.json --gml-dir /tmp/tcs-selected-gml
```

These calculations handle the explicit NoEdge constraints and the stated
valence filter. They do not perform the full global chemical-constraint
completion or construct a complete chemical CTMC or moment hierarchy.
