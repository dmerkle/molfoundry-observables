# Eight sites: one ODE across three topologies

The [input](../examples/eight-slots/model.json) declares exactly eight C sites,
each with two distinct `next` neighbours, and the
[local ownership axioms](local-slot-experiment.md). The three binding rules
preserve every vertex, ownership edge and `next` edge. Only occupancy changes.
The selected query `P0` counts occupied sites.

A finite simple two-regular site graph consists of cycles. With eight sites,
the possible topologies are `8`, `4+4` and `3+5`. Connectivity is not assumed.
Ownerless P/L vertices are permitted but inert; the supported patterns require
complete owned sites. Finite state counts below refer to the eight-site core.

## Counting identities

The [reducer](../src/mfobs/invariants/eight_slots.py) first applies the local
ownership and zero-pattern tests. Unsupported partial ownership or disconnected
site patterns remain `UNKNOWN`.

**Cycle lengths.** A cycle of length m fits when the remaining `8-m` sites
form no component or cycles of length at least three. Thus m is in
`{3,4,5,8}`. Complete six- and seven-cycle patterns count zero, as do patterns
requiring more than eight distinct C vertices.

**Path completion.** A path match containing an oriented site edge extends
uniquely around its cycle component: at each step, degree two determines the
next vertex. Its count is the sum of its possible cycle completions, with
coefficient one. A single-site pattern instead has two oriented completions,
giving coefficient `1/2`. New sites have unspecified occupancy. NoEdge tests
between original C vertices can exclude individual completions.

**Free sites.** Let H_(m,S) count oriented embeddings of a complete m-cycle
whose positions in S must be occupied; other positions are unrestricted.
For required-free positions F, use the exact Boolean identity

```text
product(i in S) x_i * product(j in F) (1-x_j)
  = sum(T subset F) (-1)^|T| product(i in S union T) x_i.
```

NoEdge remains active during MolFoundry composition; this identity rewrites
the resulting count function. Rotations and reflections identify equal
monomial counts. An internal key such as `C8:00000001` denotes one required
occupied position; `0` means unrestricted, not free.
In particular,

```text
P0 = (1/2) sum(m in {3,4,5,8}) H_(m,{0}),
1  = (1/16) sum(m in {3,4,5,8}) H_(m,empty).
```

The second identity follows because each m-cycle has `2m` oriented embeddings
and the total number of sites is eight. These are graph-count identities,
not factorizations of random moments.

## Composition and closure

The [adapter](../src/mfobs/models/eight_slots.py) supplies these identities to
the shared linear worklist. It composes both orders with all three rules,
including empty and disconnected overlaps. Seven incidence upper bounds
implied by Ω prune impossible partial overlaps: one owner/owned partner for
each relevant P/L relation and at most two `next` neighbours per C. Exact
degree two is used during cycle completion, not imposed on partial patterns.

Exact rational elimination of formal component-polynomial coefficients
selects new directions. No host occupancy or orbit vector is evaluated during
synthesis. The resulting 48-coordinate expectation ODE is closed. Subsequent
evaluation finds semantic rank 47: a triangle occurs exactly when a five-cycle
occurs, so `H_(5,empty) = (5/3) H_(3,empty)`. This additional inventory identity
was not used by synthesis; minimality is not claimed.

## Independent checks and initial values

After synthesis, the checker enumerates all 256 occupancies per topology:

| Topology | Labelled occupancies | Symmetry classes |
|---|---:|---:|
| 8 | 256 | 30 |
| 4+4 | 256 | 21 |
| 3+5 | 256 | 32 |
| Total | 768 | 83 |

It checks 36,864 pattern counts and 110,592 rate-separated generator rows
exactly. The 83 classes are not a target basis for symbolic ODE discovery.
The separate finite-quotient event construction does use their enumeration.

For a concrete initial graph G, count every discovered pattern to obtain
`mu_i(0) = count(P_i,G)`. The same coefficient matrices apply to all three
topologies; their initial vectors differ. Trajectories are checked against
independent 256-state generators. Initial graphs do not feed back into synthesis.

## Run and inspect

From the package directory:

```sh
python -m tutorials eight-sites --out /tmp/eight-sites
python scripts/run_eight_slots.py
python scripts/check_eight_slots.py
```

The first command saves a readable calculation. The other two regenerate analysis outputs and independent checks. Inspect the
`observable definitions`,
`ODE`,
`verification` and
full record produced with `--details`.
This argument relies on the supplied size-eight and ownership axioms; it is
not a closure result for arbitrary ring sizes.
