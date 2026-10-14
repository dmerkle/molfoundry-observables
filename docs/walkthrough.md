# From a graph query to an ODE

Start with the [graph tutorial](../tutorials/graph_closure_tutorial.py).
For a manifest-based example, the [gated path](reusable-family.md) supplies
one vertex query, graph rules and a recursive component language:

```sh
python run.py inspect examples/gated-path/model.json
python run.py run examples/gated-path/model.json --out out/gated-path --progress
python run.py verify examples/gated-path/model.json --result out/gated-path/result.json --conditions
python run.py ode out/gated-path/result.json
```

`run` exports the discovered coordinates and rate-separated equations.
`verify` checks the calculated generator against direct graph applications.
The [backend guide](public-api.md) explains which all-state counting and
preservation arguments are supplied by the graph-family backend.
