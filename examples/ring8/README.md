# Fixed ring: exhaustive finite analysis

This example computes observable dynamics on one fixed eight-site ring.
It uses exhaustive state evaluation to reduce counts and verify the resulting
30-coordinate ODE and 196-event model. The [Eight sites example](../eight-slots/README.md)
instead synthesizes equations symbolically for three possible core topologies.

The actual input is [`examples/ring8/model.json`](../../examples/ring8/model.json).
It lists the three conditional-DPO GML rules in `rules/` and the **sole**
selected observable `P0` in `observables/occupied.gml`. Its
`state-language.json`, `reduction.json` and `parameters.json` state the
eight-site assumptions, dihedral quotient and rates. There are 256 labelled
occupancy states and 30 dihedral classes; these are a complete finite domain,
not a length-cutoff extrapolation.

```sh
PYTHONPATH=src python -m mfobs.example_cli ring8 --out /tmp/ring8.json
```

The readable path is [`run.py`](run.py) →
[`worklist.run`](../../src/mfobs/core/worklist.py) for the linear
commutator expansion → [`ring8_common.analyze`](../../src/mfobs/models/ring8_common.py)
for the finite rule projection and generic pJC compiler. In the JSON, read
`linear` for the 30-coordinate closed ODE basis, then
`pjc.finite_language_compilation.event_classes` for the 196 derived jump
classes and `pjc.operator` for the first-order EMGF operator.
`finite_rule_projection` records the complete quotient checks.
The fixed model and symmetry reduction are supplied assumptions; the basis,
event classes and operator are computed rather than inserted into the input.

The selected identity can be replaced through the common interface, without
editing `model.json`:

```sh
PYTHONPATH=src python -m mfobs.example_cli ring8 \
  --observable P0=examples/ring8/observables/occupied.gml \
  --out /tmp/ring8-explicit-query.json
```
