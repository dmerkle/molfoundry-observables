# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Eight sites
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/eight_slots_analysis.py) · [Model inputs](../../examples/eight-slots/model.json)
#
# One occupied-slot query, three possible topologies, one symbolic ODE.

# %%
from pathlib import Path
import sys
PACKAGE = next(p for p in (Path.cwd(), *Path.cwd().parents)
               if (p / "src/mfobs").is_dir())
NOTEBOOKS = PACKAGE / "notebooks"
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / "src"), str(PACKAGE)]
import matplotlib
matplotlib.use("module://matplotlib_inline.backend_inline")
import visual_view as visual
import visual_suite as suite

run = visual.Run("visual-05_eight_slots")

# %% [markdown]
# ## 1 · Rules & state space
#
# Each site owns one binding slot. Eight sites form a simple undirected graph of degree two, giving the possible topologies 8, 4+4 and 3+5.

# %%
from mfobs import analyze as analyze_model
from mfobs.applications import EightSlots
backend = EightSlots()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
visual.rules(list(inputs.rule_gmls().items()))
suite.topology_picture((0,) * 8)

# %% [markdown]
# ## 2 · Inputs & semantics
#
# P0 counts occupied own P–L bound edges. This is the sole input observable.

# %%
query = inputs.query_gmls()["P0"]
visual.rules([("P0 · occupied slot", query)], identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# Compose the rules with the query and reduce the resulting counts using the graph-family invariants.

# %%
MAX_EXPANSIONS = 100
result = analyze_model(inputs, backend=backend, max_expansions=MAX_EXPANSIONS).record()

print("Closure:", result["status"])
print("Coordinates:", len(result["basis"]))
suite.gml_table(result["identity_patterns"])

# %% [markdown]
# ## 4 · ODE
#
# The generated system has 48 coordinates.

# %%
suite.coefficient_equations(result, names=('P0', 'P1'))

# %% [markdown]
# ## 5 · Initial states
#
# Now choose occupancy, after synthesis. The same eight-bit pattern is placed on each topology; the complete coordinate vector is evaluated separately.

# %%
BITS = (0, 0, 0, 0, 0, 0, 0, 0)
initial = suite.eight_initial(result, BITS)
suite.topology_picture(BITS)
shown = [name for name in result["basis"] if any(row[name] for row in initial.values())]
visual.table([{"Topology": tag, **{name: row[name] for name in shown}}
              for tag, row in initial.items()])

# %% [markdown]
# ## 6 · Dynamics
#
# Integrate the same generated 48-coordinate system with three different initial vectors.

# %%
T = 30.0
rates = {'basal_bind': 0.08, 'cooperative_bind': 1.2, 'unbind': 0.4}
curves = suite.eight_curves(result, initial, rates, end=T)
controls = suite.rate_controls(lambda rates, **opts: suite.eight_curves(result, initial, rates, **opts), rates, curves['options'])

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
