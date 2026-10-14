# DNA walker: MolFoundry / FRAGGER comparison

The [notebook](../../../examples/06_dna_walker.ipynb) recomputes the MolFoundry equations and compares them with the retained FRAGGER output below. [Inputs and execution](../../../../examples/dna-walker/README.md).

23 September 2026. **PASS.** The comparison uses the two-legged DNA walker
from Danos et al., *Rate Equations for Graphs*, CMSB 2020, Section 1.1.

Starting with only the extended-walker observable, the MolFoundry composition
worklist discovers the compressed-walker count and closes with two equations.
Output from the pinned FRAGGER core agrees on both raw incidence equations
and all eight individual rule/observable generator contributions, with zero
discrepancies. After the counting normalization described below, the two
configuration equations and stationary mean velocity agree symbolically with
the published result. No timing or speedup comparison is claimed.

The raw compressed pattern has two injective matches per physical compressed
configuration because its two parallel leg edges can be permuted. Writing its
raw count as `Craw`, the published configuration count is `C=Craw/2`.
The raw cross-tool comparison uses unit per-match multipliers. For the published
configuration-rate convention the extension rules use multiplier 1/2 and the
compression rules multiplier 1. The resulting equations are

```text
L E = (kFE+kBE) C - (kFC+kBC) E
L C = -(kFE+kBE) C + (kFC+kBC) E
```

Together with `C+E=1`, the stationary velocity derived from the literal
half-step rule displacements is

```text
V = (kFE*kFC-kBE*kBC)/(kFE+kBE+kFC+kBC)
```

and its symbolic difference from the expression printed by Danos et al. is zero.

The declared family is one simple directed circular DNA scaffold of arbitrary
finite length at least five, carrying one two-legged walker. The legs attach
either to the same site or to adjacent sites. Edge identity is retained through
an incidence encoding. The reducer accepts complete two-leg walker patterns
with connected DNA-path context shorter than five sites; at a fixed walker
position such context extends uniquely around every permitted cycle. Unsupported
pattern forms are refused rather than assigned an inferred count.

Independent post-synthesis checks cover circular scaffolds of lengths 5--12:
136 hosts, 408 literal multigraph rule applications, 1,088 generator-image
checks, 272 original/incidence observable-count checks and 4,536 checks of
individual projected patterns. All pass. Synthesis itself performs no host-state
enumeration and records 16 ordered composition calls and 800 raw composites.

The exact MolFoundry wheel supplied for the run has SHA-256
`c471f88069fd612f447e8402c7824f467e87fdd25b422c059434768a6715929a`;
42 installed package files were verified against it. FRAGGER used the unmodified
core at `0d9595c5d6e8918fc605f7d65778568a74f5e2d4`. A fresh GitHub Actions
execution of the original bimotor example and the comparison exports is recorded
as run 35866705340, job 107200017997. The exact result fields are in
[result.json](result.json), and [provenance.json](provenance.json) records the
runtime and source fingerprints.
