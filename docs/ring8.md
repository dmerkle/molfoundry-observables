# Fixed ring: exhaustive finite analysis

This is the eight-site binding model restricted to a single ring. It differs
from the [symbolic eight-site analysis](eight-slot-topologies.md), which uses
structural identities across three topologies without enumerating occupancies.

## Inputs and matching

Eight C/P/L site gadgets give 24 fixed vertices. Each site's optional P--L
`bound` edge records its occupancy. There are 256 labelled occupancy states
and 30 classes under rotations and reflections.

The three rules are `basal_bind`, `cooperative_bind` and `unbind`. A free site
with b occupied neighbours binds at rate `k_on + b*k_coop`; an occupied site
unbinds at rate `k_off`. Every neighbour match contributes separately.
C/P/L and `cp`/`cl`/`next`/`bound` are abstract graph labels.

## Discovery and exact check

The only selected query is the occupied-site count `P0`. MolFoundry composition
generates further graph patterns; their exact count vectors on the 30 classes
determine linear dependence. This finite-language enumeration is part of the
backend. A target observable basis is not supplied.

The retained result has 30 coordinates. Its coefficients keep rate channels
separate. The independently constructed bit-flip generators satisfy

```text
Q_r B = B A_r^T,
```

where B stores observable functions as columns over state classes and A_r
acts on the column of their expectations. The checker also verifies lumpability
on all 256 labelled states. The certificate is exhaustive for this finite
family; it makes no claim about arbitrary ring sizes.

From the package directory:

```sh
python run.py run examples/ring8/model.json --out /tmp/ring8
python run.py verify examples/ring8/model.json \
  --result /tmp/ring8/result.json --multiplicities
```

Inspect `ode.txt`, `observables/` and `coefficients-by-rate.csv` in the output.
The [retained result](../examples/ring8/README.md) records the
actual discovery steps and composition totals.
