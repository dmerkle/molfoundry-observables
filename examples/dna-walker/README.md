# DNA walker input

`input.json` contains exactly the four original directed-multigraph rules and the selected extended query from the [recovered FRAGGER export](../../notebooks/visual/evidence/dna-walker/recovered/walker-export.json). It contains no auxiliary query, basis or equations. Tests compare these input fields exactly.

Run from the package directory:

```sh
PYTHONPATH=src:. python -m tutorials.dna_walker_inspection
```

The [backend](../../src/mfobs/models/dna_walker/counting.py) declares one two-legged walker on a simple directed cycle of arbitrary length ≥5, with coincident or adjacent leg targets. It certifies DPO dangling and reduces short connected scaffold paths by unique extension; unsupported forms are not treated as zero. Rule preflight requires the scaffold and motor to be preserved. The shared Algorithm 1 loop discovers residual coordinates. Configuration names, automorphism normalization and stationary motion are computed afterwards.

The [notebook](../../notebooks/examples/06_dna_walker.ipynb) adds drawings, integration, independent literal multigraph checks and comparison against retained FRAGGER output. FRAGGER output is retained comparison data; the notebook does not rerun FRAGGER.
