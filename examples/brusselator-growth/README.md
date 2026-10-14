# Bond-isomer Brusselator with growing lactone chains

This additional example combines a molecular encoding of the Brusselator with
X-catalysed lactone growth. Its event model is exact on the declared family;
the nonlinear ODE plotted in the notebook is a mean-field approximation.

## Inputs

`BrusselatorROP().specification()` supplies Ω, rated rules and the X/Y queries.
The family contains a chain of arbitrary finite length and arbitrary finite
X/Y and monomer populations, with A/B reservoirs and discarded D/E products.
Molecules have explicit hydrogens and element/bond labels:

| Coordinate | Molecule | Formula |
|---|---|---|
| A | cyclobutane | C4H8 |
| X | 1-butene | C4H8 |
| Y | methylcyclopropane | C4H8 |
| E | isobutene | C4H8 |
| B | cyclopropane | C3H6 |
| D | propene | C3H6 |

The four atom-preserving core rules encode:

```text
A       -> X
2 X + Y -> 3 X
B + X   -> Y + D
X       -> E
```

Two more rules supply monomer and insert it into the chain, preserving the X
catalyst. These are hypothetical net reactions, not proposed elementary
hydrocarbon mechanisms. Literal inputs are in this directory; their construction
and state family are in [the molecular model](../../src/mfobs/models/brusselator_growth_chemistry.py).

## Analysis and dynamics

The X/Y queries discover higher factorial moments. For example, the X generator
contains `(N_X)_2 N_Y`; its generator introduces further moments. Reaching the
linear worklist budget therefore reports UNKNOWN rather than a closed ODE.

The event branch projects complete molecular components to population rules and
uses the shared composition compiler. A separate reaction-channel calculation
checks the result. The notebook shows the generated events and compares SSA
with a factorized deterministic approximation. For a=1, b=1.8 gives damped
oscillations and b=2.2 gives a limit cycle in that approximation. A stochastic
ensemble mean is not claimed to be a periodic orbit.

## Run

Install from the [root instructions](../../README.md#install-and-open), then run:

```sh
python -m tutorials brusselator
python -m tutorials brusselator --out /tmp/brusselator
jupyter lab notebooks/examples/10_brusselator.ipynb
python -m pytest -q tests/test_brusselator_growth.py
```

The CLI saves files only with `--out`; the notebook displays rules, equations
and dynamics in its cells. The [interactive website](https://stochmech.algochem.techfak.de/experiments/brusselator-growth/index.html)
is maintained separately.
