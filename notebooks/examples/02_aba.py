# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # ABA on a growing path
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/aba_analysis.py) · [Model inputs](../../src/mfobs/models/aba/inputs.py)
#
# Count a motif in a growing A/B path. Change `MOTIF` to try another query; `GATE=True` adds a gate that blocks appending. For the paper’s numerical experiment, use `GATE=True` and `T=20`.

# %%
from pathlib import Path
import sys
NOTEBOOKS = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'paper_view.py').is_file()))
PACKAGE = NOTEBOOKS.parent
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import paper_view as view
import visual_view as visual

# %% [markdown]
# ## 2 · Inputs & semantics

# %%
from mfobs.applications import ABA
MOTIF = 'ABA'
GATE = False
backend = ABA.motif(MOTIF, gated=GATE)
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)

# %% [markdown]
# **Algorithm 1 · Linear expectation closure**

# %%
linear = backend.linear_closure(inputs, maximum_dimension=24)
print("Linear expectation closure:", linear["status"])

# %% [markdown]
# **Algorithm 2 · Pattern support and polynomial events**

# %%
result = backend.jump_closure(inputs)
print("Polynomial jump closure:", result["pjc"]["status"])

# %% [markdown]
# ## 3 · Observable discovery · 4 · Closure status

# %%
result["linear"] = linear
support = result["support_discovery"]
print(f"{len(support['patterns'])} counts · {result['pjc']['event_class_count']} events")
print("Event model:", result["pjc"]["status"])
print("Linear expectation closure:", result["linear"]["status"])
visual.table([{"Count": name, "Meaning": p["kind"], "Word": p.get("word", ""),
               "Role": "query" if name in support["selected_queries"] else "discovered"}
              for name, p in support["patterns"].items()])
view.show_equations(result, names=support["selected_queries"])

# %% [markdown]
# The event model closes; linear expectation closure remains unresolved at the budget of 24 functions.
#
# Finite-state projection (FSP) approximates the forward equations. SSA samples the same events. The plot reports truncation bounds and sampling uncertainty.

# %% [markdown]
# ## 5 · Initial states

# %%
WORD = "ABAB"
POOLS = {"A": 2, "B": 1}
T = 4.0
RUNS = 100
rates = {"beta_A": 2., "beta_B": 1.5, "k_A": .9, "k_B": .01}
if GATE:
    rates.update(a=.6, b=.4)
initial = view.word_coordinates(result, WORD, POOLS, opened=True)

# %% [markdown]
# ## 6 · Dynamics

# %%
visual.close_controls(globals().get("controls"))
curves = visual.show_curve(
    visual.aba_curve(result, initial, rates, end=T, paths=RUNS)
)
controls = visual.aba_sliders(result, initial, rates, curves)

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
