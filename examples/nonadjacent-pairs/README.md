# Counting and condition test inputs

An atom-preserving cycloaddition diagnostic for the generic composition and
condition tests. `guarded.gml` excludes directly linked alkene pairs;
`ablated.gml` removes only those two context guards. Neither input is a
validated kinetic model.

The raw rate factor is `1/8`: eight injective atom maps represent each unordered
ethylene pair. On mixtures of ethylene and cyclobutane the propensity is
`k*n*(n-1)/2`, so the expectation hierarchy does not close uniformly over all
initial inventories.

```sh
python -m pytest -q tests/test_nonadjacent_chemistry.py tests/test_pair_depletion.py tests/test_noedge_zero.py
```

These files are test inputs, not an additional paper example. Independent
finite checks live in the tests; the synthesis uses the shared worklist.
