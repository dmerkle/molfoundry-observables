# MØD 1.0.1.100 / MolFoundry Aldol comparison

Executed on 25 September 2026 with MØD 1.0.1.100. The input script is
`calc-original.py`; `report.json` records its hash. `audit.py` wraps its result
handler to export successful rules without changing the calculation or rule
allocation order. The exported GMLs and summary record retain the results. Reference
rules were loaded after all four runs, preserving their rule IDs.

| Input and overlap implementation | After | Before | Distinct rule classes |
|---|---:|---:|---:|
| Plain, Boost=True | 9 | 18 | 21 |
| Plain, Boost=False | 13 | 18 | 22 |
| `:` absence encoding, Boost=True | 8 | 16 | 19 |
| `:` absence encoding, Boost=False | 10 | 18 | 23 |

After means reaction then identity; before means identity then reaction. In the
original output these are respectively count21 and count12. The Boost flag,
not the version number alone, selects the restricted versus complete enumeration.

## Plain comparison

`plain-complete/r978.gml` is isomorphic to MolFoundry class
`new-f7c8786fa6b1`, verified using MØD's Rule.isomorphism. Every one of the 31
accepted contributions agrees by rule class and per-order multiplicity.
The 13-versus-12 difference is against the published TCS table after its manual
alias names have been merged. MØD 1.0.1.100 and MolFoundry agree.

## Absence-condition comparison

The second model uses ordinary `:` edges as absence sentinels with valence
weight zero. MØD sees graph edges in these GML files, not native NoEdge conditions.
`encoded-noedge-complete/` contains these original exported products.

For each product, the explicit translation removes `:` edges and converts each
left/context sentinel pair into an input `constrainNoEdge`. No right-only sentinel
edges occur in these products. The translated rules are in `translated-noedge/`.
All **23 conditional classes** and all **28 contributions** match fresh native
MolFoundry conditional composition, with after10/before18 and 18 nonzero signed
terms. The native products and full fresh record are in `molfoundry/`.

`cross-tool-comparison.json` gives every original/translated/native GML path,
canonical identity, and multiplicity check. The translation was checked for this
example's resulting rules; it is not a general theorem about complement encodings
or evidence that MØD natively supports arbitrary NACs.

The recorded scripts contain absolute paths from the execution environment.
The notebook recomputes MolFoundry compositions and checks them against these
saved MØD exports; it does not rerun MØD.
