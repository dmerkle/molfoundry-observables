# Exact linear compression

Compression accepts a closed exact result and computes an invariant subspace
containing its selected queries. It preserves the source result. From the
package directory:

```sh
python run.py run examples/gated-path/model.json --out /tmp/gated-path
python run.py compress /tmp/gated-path/result.json --out /tmp/gated-path-compressed
python run.py ode /tmp/gated-path-compressed/result.json
```

The gated-path result already has three necessary coordinates, so its dimension
stays three. Compression can remove unnecessary coordinates from other closed
inputs. It cannot complete an unresolved discovery run.

For source coordinates `f` and rate-separated generators `A_r`, an exact
rational worklist closes the selected rows under `v -> v A_r`. Channels sharing
a rate name are added first; distinct rate names remain independent. If `B`
contains the resulting independent rows, the compressed coordinates are `g=B f`.
The implementation checks

```text
B A_r = C_r B,
seed row = seed_reconstruction row times B.
```

The output includes the transformed ODE, `transformation.csv`,
`basis-combinations.json` and `seed-reconstruction.json`. Compressed coordinates
are linear combinations of source counts; they need not be single graph
patterns. Available initial values are transformed by the same matrix.

This is the smallest common invariant rational subspace containing the selected
queries, conditional on the source coordinates being independent functions on
Ω. Source closure and preservation of Ω remain assumptions. Numerical rate
specialization or further count identities can permit additional reduction.

[Implementation](../src/mfobs/invariants/compress.py) ·
[Tests](../tests/test_linear_compression.py).
