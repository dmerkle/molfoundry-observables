# Rémy trees

The example reproduces §5 of Nicolas Behr, *On Stochastic Rewriting and
Combinatorics via Rule-Algebraic Methods* (2021). Published equations are used
in checks, not supplied to discovery.

## Encoding and Ω

Ω consists of valid finite planar rooted binary trees. MolFoundry GML edges are
undirected, so each directed edge becomes an auxiliary node labelled `I`, `L`
or `R`, connected to its source by `src` and its target by `tgt`. Object vertices
have label `V`. This retains direction, edge type and match multiplicity.

The six insertion rules select an I/L/R edge and place its subtree on the left
or right of a new branch. They delete the selected edge-node and its two
incidences, preserve its endpoints, and create the branch, leaf and ordered
child edges. No NoEdge guard is needed on valid trees.

Composition rejects nonfunctional incidence ports, cycles and incompatible
edge configurations. Positive child-totality and connectedness belong to Ω;
they are used by the tree counting backend.

## DPO projection and count reduction

The source model uses SqPO. On Ω, the selected edge-node has exactly the two
matched incidences, so DPO has the same input matches and targets.

MolFoundry's DPO jump projection is `I ← K → I`. The tree backend first computes
it, checks that only complete edge-nodes are deleted, then reduces its count to
the full-input identity `I ← I → I`. Their match counts agree on Ω: there are no
unmatched incidences and recreation returns an isomorphic host. This is a count
identity on valid trees, not an equality of arbitrary rule spans.

The reducer rejects incomplete deleted edge-nodes, object deletion and
unsupported conditions. A separate regression covers 384 generator/composite
rules, including canceled terms.

## Discovered counts and events

`P_k` counts internal vertices with a left spine of at least k internal vertices.
Only `P3` is selected. Composition and exact linear-span reduction discover

```text
Q P3 = -4 P3 + 4 P2,
Q P2 = -2 P2 + 3 P1,
Q P1 = 2 E,
Q E  = 4 E.
```

`E` counts edges. Names are assigned after structural recognition of the
composed patterns. These identities are derived as canonical rule expressions,
not fitted on enumerated trees.

For an insertion edge, let d be its subtree's left-spine depth and a its number
of consecutive left-child ancestors. The count of edges with `d>=r, a>=s` is
`P_(r+s)`, by a tree bijection. Exact finite differences partition these threshold
counts into `(d,a)` bins. Grouping their observable jumps gives four event
classes; higher-spine terms cancel, leaving activities in `E,P1,P2,P3`.

The shared worklists perform observable-span and operator-span closure. The
backend supplies the incidence-tree constraints and the bijection above.
The unpartitioned source rules have constant jumps only in `E,P1`; the deeper
spine counts require the event partition.

Independent enumeration checks 2,056 trees and 66,196 rewrite events through
size eight. This finite check supplements the all-tree counting argument.

## EMGF and distributions

The derived EMGF operator agrees with Eq. 57g. Its jump monomials differ from
printed Eq. 58; the comparison records this discrepancy. Integer propagation of
the jump operator yields the joint `(P2,P3)` distributions without enumerating
all trees, checking total weight `(2n)!/n!` before normalization.

From the package directory:

```sh
python -m tutorials remy
python scripts/build_remy_trees.py --out /tmp/remy-trees
pytest -q tests/test_remy_trees.py
```

`--skip-raw-compositions` omits the raw expression exports when rebuilding the
report. [Inputs](../src/mfobs/models/remy/inputs.py) ·
[Counting](../src/mfobs/models/remy/counting.py) ·
[Events](../src/mfobs/models/remy/events.py).
