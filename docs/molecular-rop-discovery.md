# Direct molecular ROP observable discovery

The [molecular ROP analysis](../tutorials/rop_analysis.py) uses native
atom-resolved rules and query patterns. Change `MOTIF` or `GATED` in that
file, or use the [notebook](../notebooks/examples/07_molecular_rop.ipynb).

```sh
python -m tutorials rop
python -m tutorials rop --out /tmp/rop --details
```

Algorithm 1 composes rules and reduces their count images directly.
Algorithm 2 independently discovers support and extracts events. Both use
[MolecularReducer](../src/mfobs/models/molecular_rop/counting.py);
typed-path constructors do not participate in molecular discovery.

## Actual computation and retained evidence

1. Load atom-resolved feed, alcoholysis and optional Pd-N gate rules. Only the
   selected chemical motif is a query. Feed markers denote boundary operations.
2. Compose both orders with common overlaps, including empty/disconnected ones;
   transport conditions, retain multiplicities, cancel signed terms and apply
   DPO jump projection with the interface retained.
3. Reduce the resulting **molecular** components. Recognise free lactones,
   whole incorporated repeats, OH-anchored endings, the unique OH and gate.
   Sequence names are assigned at this point, after composition. New factors
   enter the shared support worklist; no suffix inventory is prepared.
4. Check source identity and the model-specific uniform-match argument. Partition
   Boolean ending/gate counts to obtain constant joint jumps and polynomial
   activities. Use the existing polynomial generator, expectation and EMGF code.
5. Validate against independent molecular rule application and, separately,
   the typed-path derivation. Neither reference supplies target equations.

With `--details`, the saved support trace records the operands, signed
composites, conditions, count reductions and discovery origins.

## Supplied all-chain arguments (not inferred from tests)

The state family is exactly the existing constructors: one benzyl-initiated,
unbranched chain, finitely many free valerolactone/caprolactone molecules,
feed markers, and optionally one asymmetric Pd-N complex. Hydrogens are explicit.
No branching, backbiting, transesterification or other chemistry is claimed.

**Local overlap pruning.** The empty and one-unit chain constructors, free
monomers and both gate states give all possible atom neighbourhood signatures
(labelled incident bonds). Interior ester oxygens have the same two single C
ports as the ester oxygen in a one-unit chain; all remaining internal atom
neighbourhoods also occur there. Thus these signatures cover every chain length
and monomer sequence. The source rules preserve this family. These profiles are passed to MolFoundry as public incidence constraints.
The kernel uses their necessary local bounds to prune overlaps. Nonlocal
restrictions are checked by public predicates; no private engine code is used.

**Zero connected patterns.** A connected pattern with n vertices that embeds
in a polymer can meet at most floor((n-1)/6)+2 consecutive repeat blocks. Between
non-boundary blocks, every repeat contributes at least six backbone edges;
a connected injective pattern spanning them must contain those edges. Removing
unmatched prefixes/suffixes therefore places its image in a chain of at most
that many repeats. Enumerating all words up to this derived bound covers every
root, interior, end and short-chain case. Free molecules, markers and both gate
states complete the host-component cover. Absence of a *positive* embedding in
this cover proves count zero on all constructor states, also when additional
negative conditions are present. Finding an embedding does not certify a
formula: an unrecognised positive shape returns UNKNOWN. The depth budget is 8.
This finite local-context argument is distinct from bounded transition testing.

**Recognised counts and conditions.** Carbonyl orientation and the spacer
length distinguish the A/B units. A whole repeat/end pattern has one injective
map per occurrence, including overlapping occurrences; a free cyclic ester has
one oriented map per free molecule. OH occurs only at the unique chain end.
The recognised organic patterns are induced on their images in the constructor
family, so their additional internal absent-edge predicates are redundant unless
contradictory. The Pd-N absence condition is retained as the Boolean gate count.
It is not erased as an organic-pattern identity.

**Factorisation and deletion.** Only one connected chain factor and at most
one free-pool factor may occur in an accepted product. They necessarily occupy
different host components, as do the gate and distinct feed markers; cross-
component absent-edge tests are then automatic. Multiple chain/pool factors
are refused rather than assuming independence or losing falling-factorial
corrections. Two copies of the same unique feed marker have zero injective
count. Molecular rules preserve existing atoms; a projected rule containing
vertex deletion is explicitly unsupported. Thus no DPO dangling condition is
silently discarded by converting its input into an identity count.

**Uniform jumps.** At a fixed state each insertion match selects the unique
OH and one oriented free lactone. Different copies give isomorphic molecular
successors but remain distinct rate contributions. Feed and gate matches are
unique. Consequently each source's generator image divided by its activity
is its coordinate increment on enabled states. This is an additional source-
family argument, not a consequence of first-moment closure. Source identities
are checked before the shared finite Boolean partition algebra is used.

## Validation and limitations

Separate checks in [checks/molecular_rop.py](../src/mfobs/checks/molecular_rop.py)
compare activities, generator images and joint jump/rate distributions with
direct MolFoundry applications and independent concrete counts. They also
compare with the typed-path model under the explicit graph-count correspondence.
These checks run after synthesis; they do not supply equations or prove the
all-chain arguments above.

Support closure, finite events and linear expectation closure have separate
statuses. Budget failure retains the pending queue and exact operands in a
minimal reproducer, without exporting a partial event operator. The candidate
budget is per composition order, after sound profile pruning:
`--molecular-overlap-limit` (default 2,000,000). Support pattern/expansion and
linear-dimension budgets retain their existing options. A wall-clock limit,
if needed, is enforced by the reproduction driver in a separate process.

The reducer supports this constructor family and its stated condition
fragment. Unsupported positive shapes, deletion interfaces or conditions
remain unresolved.
