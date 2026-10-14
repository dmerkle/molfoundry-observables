# Counts on a chain and finite components

The `chain_mixture` adapter reads a recursive chain language, two finite
controller states and a repeated finite species. It discovers coordinates
from composed graph patterns, without a supplied target basis or equations.
The [gated path](reusable-family.md) uses this adapter.

For a connected pattern with k vertices, a match occupies at most k consecutive
chain cells. The counter enumerates these finite footprints, distinguishing the
left boundary, right boundary, interior, both boundaries and the exceptional
zero-length chain. Interior footprints translate along the chain; short-length
corrections are retained exactly. No finite-length polynomial fit is used.

Disconnected patterns require joint injectivity. The counter corrects products
for vertex identifications and cross-component NoEdge conditions; see
[joint chain counting](joint-chain-counting.md). Components assigned to repeated
finite molecules contribute falling factorials in the number of copies.
The two controller states are combined in the Boolean coordinate x, with x²=x.

These identities depend on the declared grammar and supported conditions.
Rule preservation of that language is a separate backend obligation. Unsupported
conditions and exhausted budgets return `UNKNOWN`. Concrete graph checks test
the implementation; they do not prove the all-size language contract.

The implementation is in [chain_mixture.py](../src/mfobs/models/chain_mixture.py),
[chain.py](../src/mfobs/invariants/chain.py) and
[counting.py](../src/mfobs/invariants/counting.py).
