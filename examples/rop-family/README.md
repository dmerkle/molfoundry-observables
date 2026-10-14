# Lactone ring-opening polymerization

A growing chain records A/B insertions. A is δ-valerolactone and B is
ε-caprolactone. Each insertion consumes one free monomer; constant-rate feeds
supply the pools. An optional missing-edge gate controls insertion.

The molecular insertion is an atom-preserving net alcoholysis model. The
Pd/P/N gate is a software model of an active relational condition, with
illustrative rates. It is not attributed to an experimentally validated
catalyst mechanism. The chemical modeling reference is
[Neilson and Bielawski (2013)](https://doi.org/10.1039/C3CC42424C).

## Select a query

From the package directory:

```sh
python run.py rop-analysis --motif ABA --linear-limit 24 --out /tmp/aba.json
python run.py rop-analysis --motif ABBA --gate --out /tmp/abba-gated.json
```

The first command selects only the ABA count. Compositions and counting
reduction discover the pools and terminal patterns `M_A`, `M_B`, `S_AB` and
`S_A`. A unique path end and interchangeable free tokens justify uniform
match jumps. This family argument is necessary for event extraction in
addition to polynomial support closure.

The typed path analysis and atom-resolved molecular composition have separate
inputs and provenance. Their correspondence and independent host checks are
available in the [molecular notebook](../../notebooks/examples/07_molecular_rop.ipynb).
[Support discovery](../../docs/rop-support-discovery.md) documents the result
fields, counting identities and budget exits.

The [ABA notebook](../../notebooks/examples/02_aba.ipynb) accepts motif
changes and computes numerical dynamics from the generated events.
[Interactive website](https://stochmech.algochem.techfak.de/experiments/rop-family/index.html)
to inspect rules, highlighted matches, expansions and curves.
