# Brusselator--ROP: same query interface, two distinct conclusions

The default selected GML identities are `X` and `Y` in
[`examples/brusselator-growth/observables/`](../../examples/brusselator-growth/observables/).
The six literal rules are in the neighbouring `rules/` directory.

```sh
PYTHONPATH=src python -m mfobs.example_cli brusselator --out /tmp/brusselator.json
```

The common linear worklist reaches its dimension limit, so it does not claim a
finite expectation ODE.  The composition-derived four-coordinate event model
does have polynomial jump closure and therefore supplies the finite jump table
and EMGF operator.  This is the same separation of Algorithm 1 and Algorithm 2
used by the other examples.
