# Local site identities

This backend uses the binding rules without fixing the number of sites or
their topology. The [eight-site reducer](eight-slot-topologies.md) adds size
eight and degree two; the [fixed-ring backend](ring8.md) enumerates one ring.

## Supplied graph family

Finite simple undirected graphs satisfy these ownership axioms:

- Each C owns exactly one P through `cp` and one L through `cl`.
- Each P/L has at most one C owner and one `bound` partner.
- A `bound` edge joins only the P and L of the same owner.
- Other edges are `next` links between C vertices.
- The degree-bounded variant permits at most two `next` neighbours per C.

Ownerless P/L vertices are inert. The axioms are supplied, not discovered.
The inputs are `examples/local-slots/model-local.json`, `model-degree2.json`
and `model-no-coop.json`; the last removes cooperative binding.
`max_pattern_sites` is a reduction budget, not a restriction on host size.

## Exact structural reduction

The [local normalizer](../src/mfobs/invariants/local_slots.py) proves zero
counts for conflicting ownership, multiple bindings, bindings between distinct
owners and an edge contradicted by NoEdge. It removes a NoEdge test only when
types or ownership already exclude that edge. NoEdge on a site's own P/L pair
remains essential: absence of a positive `bound` edge does not by itself mean
the site must be free. Unsupported partial ownership remains `UNKNOWN`.

For a complete owned-site pattern, each injective map of its C vertices
determines exactly one map of their P/L vertices. The count factor is one;
automorphisms retain their multiplicity. A symmetric `00` or `11` edge
pattern therefore counts both orientations, while `01` counts one per mixed
edge. Accepted jump spans preserve their input vertices, so converting them
to identity counts introduces no dangling-condition omission.

The preservation check reads the effective MolFoundry rules. Binding adds the
own P/L edge under NoEdge; unbinding deletes it. Vertices, ownership and `next`
edges are preserved, including the optional degree bound. MolFoundry may add
the created-edge NoEdge test itself; its absence from handwritten GML alone
does not demonstrate an unguarded effective rule.

## Discovery

The [adapter](../src/mfobs/models/local_slots.py) uses the shared worklist and
all-overlap MolFoundry composition. Starting from occupied sites `P0`, it
discovers free sites `P1` and mixed adjacent sites `P2`:

```text
d E[P0]/dt = k_on E[P1] + k_coop E[P2] - k_off E[P0].
```

Further compositions introduce longer paths, triangles and, without a degree
bound, branches. The bounded cooperative searches remain `UNKNOWN`.
Removing cooperation closes with the occupied/free pair, for every graph in
the declared local family:

```text
d E[B]/dt =  k_on E[F] - k_off E[B],
d E[F]/dt = -k_on E[F] + k_off E[B].
```

## Separate nonclosure argument

A budget stop is not a nonclosure proof. A separate argument applies to the
union of all cycle sizes when the cooperative rate c is positive.
Write `N=sum_i x_i` and `m_S=product(i in S) x_i`. The bit-flip generator is

```text
Q m_S = -k_off |S| m_S
        + k_on sum(i in S) (m_(S\{i}) - m_S)
        + c sum(i in S, j adjacent to i)
            (m_((S\{i}) union {j}) - m_(S union {j})).
```

Only the negative cooperative terms raise Boolean polynomial degree. For
`n >= max(3,r+2)`, the degree-(r+1) part of `Q^r N` is

```text
(-2c)^r sum(S a proper connected cycle interval, |S|=r+1) m_S.
```

Each longer interval has exactly two predecessor intervals, obtained by
removing either endpoint; lower-degree terms cannot reach the new highest
degree. Thus `N,QN,...,Q^R N` are linearly independent on sufficiently large
cycles, for every R. No common finite-dimensional generator-invariant linear
space containing N exists over all ring sizes. This does not exclude closure
on a fixed ring, approximations or special initial distributions.

## Reproduction

```sh
python scripts/run_local_slots.py --model all --max-expansions 7
python scripts/check_local_slots.py
python -m pytest -q tests/test_local_slot_*.py
```

The independent checker runs after synthesis and compares graph counts and
generator rows. Its [reference implementation](../src/mfobs/checks/local_slots_reference.py)
also checks the Boolean-polynomial argument. Finite tests support the
implementation; the structural and degree arguments above establish its scope.
