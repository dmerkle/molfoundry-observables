# Gated path growth

This declarative example derives the expectation equations of a
two-state-modulated Yule process from labelled graph rules and one query.

## Inputs, invariant, and match normalization

The state language is the disjoint union of:

* one path on `n+1` vertices labelled `A`, with `n >= 0` edges;
* one controller with distinct labels `U`, `W`, `V`, edges `U-W-V`, and an
  optional edge `U-V`; let `x=1` mean this optional edge is absent;
* any number `M >= 0` of isolated `T` tokens. These spectators are untouched.

The single-node-tile path, both controller graphs, and repeated token are
declared in `examples/gated-path/state-language.json`. No graph-construction
function or target equations are imported by synthesis. The `chain_mixture` adapter uses the declared component language for exact
count reduction.

There are three rules:

1. `subdivide`: replace any `A-A` edge by `A-A-A`, preserving both old endpoints,
   while matching `U-W-V` and requiring **NoEdge(U,V)**. An edge has two oriented
   injective matches, so rate per match is `lambda_g/2`.
2. `close`: add the controller edge `U-V`, requiring NoEdge(U,V), rate `a`.
3. `open`: delete the controller edge `U-V`, rate `b`.

All node labels are ordinary labels. Distinct labels on the controller eliminate
controller automorphisms; no quotient of match multiplicities is silently taken.
The overlap policy is `connected=false, include_empty=true`. The forbidden
signature has degree at most two, admits only `A-A` and the three controller
edge types, and isolates `T`.

**Preservation for all sizes.** Subdividing an edge of a finite path produces
another path, adds exactly one vertex/edge, and preserves endpoint degrees.
Controller rules change only its optional edge. No rule joins distinct
components or changes a token. Each rule therefore preserves the declared
language for every `n,M`; this is a mathematical argument, not an extrapolation
from the finite checks. The `n=0` isolated path vertex has no subdivision match.
Starting with `n>=1` and positive growth/opening rates gives infinitely many
reachable path lengths and connected component sizes.

**Nonexplosion.** Path growth has total rate `lambda_g*n*x <= lambda_g*n`.
It is dominated by a linear Yule birth process, which is nonexplosive and has
finite moments of every fixed order at finite times from finite initial size.
The controller has finitely many states. For the feedback variant below its
rate is still bounded on every finite-size time interval; pathwise Yule
domination excludes accumulation of infinitely many events in finite time.

## Exact generator and discovered closure

For a function `f(n,x)`, the coordinate generator is

```
Qf = lambda_g*n*x [f(n+1,x)-f(n,x)]
     + a*x [f(n,0)-f(n,x)]
     + b*(1-x) [f(n,1)-f(n,x)].
```

This oracle is derived independently from the graph-match counts: `2*n*x`,
`x`, and `1-x`. It is evaluation code and does not guide synthesis.

Starting with only the identity on one `A`, the graph count is `V=n+1`.
Discovery gives `1 -> 2 -> 3 -> 3` with representatives `V`, `2*n*x`, and `1`.
Writing `v=E[n+1]`, `z=E[n*x]`, and `c=E[1]` gives

```
v' = lambda_g*z
z' = b*(v-c) + (lambda_g-a-b)*z
c' = 0.
```

These are exact expectation equations, not `E[n*x]=E[n]E[x]` factorization.
The third coordinate is constant; there are two nonconstant affine coordinates.
Removing only the NoEdge on `subdivide` changes its propensity to `lambda_g*n`
and the synthesized dimension to two: `v'=lambda_g*(v-c), c'=0`. Thus this
NoEdge is **dynamically active on the invariant language**, not an ambient
simple-graph redundancy. The separate NoEdge on edge-creating `close` is not
advertised as an independent novel mechanism.

## Higher-order queries

`model-degree{k}.json` selects `k` isolated A vertices, counting ordered distinct
tuples `(n+1)_k`. The rules and counting backend are unchanged. For `k>=2`, write
`R_j(n)=n(n+1)...(n+j-1)`. Then

```text
Q R_j = j*lambda_g*R_j*x,
Q(R_j*x) = b*R_j + (j*lambda_g-a-b)*R_j*x.
```

The selected query expands in `R_2,...,R_k`, giving dimension `2(k-1)` for generic
rates. The reproduction script checks degrees 2, 4, 10 and 11 against this
independent calculation. Its `--quick` option runs only degrees 2 and 4.
[Linear compression](linear-compression.md) checks invariant-subspace and query
reconstruction identities after discovery; the current worklist already reduces
linear combinations during discovery.

## A nonexplosive negative case with a proof

`model-feedback.json` changes only `close`: its context additionally contains
one `A` vertex. There are now `(n+1)*x` closing matches, so its rate is
`a*(n+1)*x`. Growth remains linear and the process is nonexplosive as above.

For `lambda_g>0` and `a>0`, no finite-dimensional linear space of functions
containing `V=n+1` is invariant under the full generator. Proof: `QV=lambda_g*n*x`.
On a degree-`j` polynomial times `x`, growth preserves the degree and closing
adds a leading `-a*n^(j+1)*x`; opening preserves degree. Inductively,

```
Q^k V = lambda_g*(-a)^(k-1)*n^k*x + terms of degree at most k-1,   k>=1.
```

The nonzero leading degrees grow strictly, so these functions are linearly
independent on `n>=0,x in {0,1}`. This also rules out finite closure by arbitrary
linear-combination observables, not merely by the current pattern heuristic.
The proof concerns functions on the declared all-size language. Degenerate
rates or special initial distributions, for example the absorbing `n=0`, can
give simpler trajectories and are not counterexamples to the claim.

The software run is bounded to five expansions. It reports `UNKNOWN`, eight
discovered functions, pending work, and **no ODE matrix**. It does not itself
prove nonclosure. Its computed partial equations are independently checked;
the theorem above supplies the nonclosure argument. No timeout is re-labelled
as a decision procedure.

## Verification and reproduction

From the project directory in the MolFoundry environment:

```
python scripts/validation/reproduce.py
python scripts/validation/reproduce.py --out out/reuse
python -m pytest -q tests/test_reusable_family.py
```

`--quick` omits degrees 10 and 11. Relative `--out` paths are resolved from the
project root; without it, the output root is `out/gated-path-checks`. Reproduction
writes only the named generated case directories, sibling `*-compressed`
directories for closed cases, and their summary. It checks compressed dimensions
and source-result hashes; `--summarize-only` requires existing matching source,
check, and compression files and never recomputes a missing result. Each case stores its input hashes, actual
rules/patterns, trace, rate-separated equations if closed, and independent
check report. The first-moment case also checks the raw composition
multiplicities independently.

The independent verifier constructs complete graphs directly, applies the
literal rules with MolFoundry matching and DPO application, recognizes
the resulting path/controller/spectator structure, counts every discovered
pattern, and checks rate-separated equations. Separately it verifies all
generated polynomial identities symbolically in arbitrary `n` and Boolean `x`.
It does not import the supplied chain grammar or the exact-counting Reducer.

The high-degree full-basis direct-host grid is deliberately small (`n=0..3`),
because brute force counting of interchangeable nodes is factorial. Those finite checks
mostly exercise lower-order discovered representatives and do **not** validate
the high-order seed at large host sizes. Separate seed-only checks construct
complete larger hosts with `k`, `k+1`, and `k+3` actual path vertices, compute
the nonzero seed count using `math.perm(actual_A_label_count,k)`, apply every
literal rule match, and check the seed generator equation. These are elementary
ordered-tuple combinatorics, not exhaustive `k!`-sized match enumeration, and
they do not claim to count every discovered high-order basis pattern. The
all-size evidence remains the exact structural counting argument and independent
symbolic identities.
