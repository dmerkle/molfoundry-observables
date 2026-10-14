# The balanced NAC example can close after a covalent tether is added

A covalent tether makes the cycloaddition intramolecular. One graph query
closes into two independent expectation coordinates.

## The modification

The rule reads two alkene bonds connected by a covalent tether and preserves
that tether during the four-carbon [2+2] bond changes:

```text
Positive reactant:  1=0-4-5-2=3    (six C atoms)
Bond changes:      0=1 and 2=3 become single;
                   create single bonds 0-3 and 1-2.
Unchanged tether:  0-4-5-2
Product safeguards: NoEdge(0,3), NoEdge(1,2)
Context conditions: NoEdge(0,2), NoEdge(1,3)
Rate per injective atom map: k/2.
```

The connected positive pattern cannot pick atoms from two separate molecules.
Each eligible molecule has exactly two matches. Therefore its population
propensity is `k*A`, not a quadratic function of the whole mixture size.
All atoms and per-atom bond-order sums are preserved; full host hydrogens are
explicit. Only C/H vertices and `-`/`=` bonds occur. There is no inlet or token
creation. NoEdge remains a predicate, not an artificial bond.

## Supplied molecular language

Arbitrary independent nonnegative copy numbers of these four finite components:

| Name | Interpretation | Formula | Explicit atoms | SMILES |
|---|---|---|---:|---|
| A | 1,5-Hexadiene precursor | C6H10 | 16 | `C=CCCC=C` |
| P | Bicyclo[2.1.1]hexane product | C6H10 | 16 | `C1CC2CC1C2` |
| B | Context-blocked diene | C6H8 | 14 | `C=C1CCC1=C` |
| D | Product of the ablated B channel | C6H8 | 14 | `C1CC23CC12C3` |

D is included so the guard-ablation comparison has the same closed language.
The species names are metadata, not atom labels or graph-generation commands.
The finite component library is supplied explicitly; it is not discovered.
There is no population cutoff. Each fixed finite closed starting inventory
still has a finite reachable state set.

The A -> P transformation is experimentally established as a selected crossed
photocycloaddition channel: Srinivasan and Sonntag, JACS 1967,
[DOI 10.1021/ja00978a042](https://research.ibm.com/publications/bicyclo-211-hexane-preparation-and-photochlorination).
The [related primary study](https://research.ibm.com/publications/mercury-lesssupgreater3lesssupgreaterplessinfgreater1lessinfgreater-photosensitized-internal-cycloaddition-reactions-in-14-15-and-16-dienes)
also reports competing isomers and dimerization. This input is not a complete
or calibrated photochemical mechanism. Light and sensitization are absorbed
into an illustrative effective rate. The hard exclusion of B, and the ablated
B -> D channel, are diagnostic modeling choices, not established chemistry.
Valence validity alone does not prove stability or synthetic accessibility.

## Synthesized equations and active condition

Only the ordinary C=C identity is an initial observable. Call its raw embedding
count `R` here (its machine ID is `B`, distinct from component B). The pipeline
finds `R=4(A+B)` and a second actual conditional jump pattern `J=2A`, with

```text
R' = -2*k*J
J' = -k*J.
```

These are equations for expectations of graph counts, with no mean-field
factorization. Equivalently, if `u` is the expected total number of alkene bonds
and `a=E[A]`, then

```text
u' = -2*k*a
a' = -k*a.
```

There is one conserved combination, `u-2*a=2*E[B]`. Thus two independent
observable functions are needed over the whole language, but there is only
one exponential relaxation mode. If the initial blocked population is known,
one can equivalently write a single affine ODE with that conserved value as a
parameter. This is not being advertised as two independent chemical modes.

This exact two-dimensional closure holds for arbitrary initial mixtures/laws
with finite relevant moments, not only one chosen batch size. The supplied
language makes connected counts linear in finitely many component populations;
this is intentionally a much simpler setting than unbounded polymer topology.
In fact, the supplied four-species language and component-local rules already
guarantee a first-moment population system of at most four coordinates. Finding
two from the chosen graph seed is a smaller exact observable closure inside
that finite ambient space, not evidence that finite closure was obtained
without a strong structural input assumption.

The context condition is active: the positive six-C pattern matches A and B
twice, but the guarded rule has two matches in A and zero in B. Removing only
the context conditions enables B -> D, giving the one-coordinate equation
`u'=-k*u`. Thus NAC ablation changes both dynamics and the required observable
basis. For example, starting with two A and two B molecules, guarded dynamics
leave four alkene bonds in the long-time limit, whereas the ablated model
consumes all eight initial alkene bonds.

## What was actually computed

- Guarded: 8,036 raw compositions, four ordered composition calls;
  one seed -> two functions -> closed.
- Ablated: 62 raw compositions, two ordered calls; one seed -> one function ->
  closed.
- Sixty independent full-molecular generator identities over twenty mixtures,
  plus direct MolFoundry species-derivation tests with finite-inventory checks.
- No target basis or ODE supplied to synthesis. Expected physical equations
  appear only in post-synthesis verification and regression tests.

The `mfobs.models.finite_components:Adapter` counts each connected
pattern on every supplied component and multiplies by that component's symbolic
copy number. It uses exact rational span tests on actual computed jump patterns.
An independent exhaustive check of local successors proves preservation for
all inventories, since connected reaction inputs lie in one component.
Disconnected noncancelled jump patterns and exhausted counting budgets give
UNKNOWN. This is not a general multivariate mixture counter or general closure
decision procedure.

## Important limitation: the NAC is not hard to eliminate

On this finite language, adding one explicit preserved `C0-H` context and
removing all NoEdges selects exactly the same maps and products. The unblocked
reactive carbon has that H; the blocked carbon does not. This positive
translation has a tested match bijection, including normalization.

Therefore this answers **“can a small balanced NAC model close?” with yes**.
It does not answer **“is this chemistry beyond FRAGGER?” with yes**. The change
that repairs closure is intramolecular kinetics, not a special NAC syntax.

## Reproduce

Use the environment from the [package setup instructions](../../README.md#install-and-open) and run these commands from the package directory.

```sh
python scripts/tethered_audit.py
python -m pytest -q tests/test_tethered_chemistry.py tests/test_tethered_synthesis.py tests/test_finite_components.py
```

The script regenerates only these dedicated example/report directories.
Results are in `out/tethered-nac/summary.json`, `guarded-result.json`,
`ablated-result.json`, the corresponding verification JSON and ODE text files.
