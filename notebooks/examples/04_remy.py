# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Rémy trees
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/remy_analysis.py) · [Model inputs](../../src/mfobs/models/remy/inputs.py)
#
# Starting from the left-spine count P3, six insertion rules discover P2, P1 and E and a closed expectation ODE. All rules share the rate κ.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
NOTEBOOKS = PACKAGE / 'notebooks'
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import visual_view as visual
from IPython.display import display, Math
import sympy as sp

# %%
from mfobs import analyze
from mfobs.applications import RemyTrees
backend = RemyTrees(order=3)
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
result = analyze(inputs, linear_limit=16, record_patterns=True, backend=backend).record()

# %% [markdown]
# ## 3 · Observable discovery · 4 · ODE

# %%
print("Discovered basis:", result["linear"]["basis"])
print("Linear closure:", result["linear"]["status"])
for step in result["linear"]["steps"]:
    rhs = sp.Add(*(sp.Rational(c) * sp.Symbol(name) for name, c in step["coordinates"].items()))
    display(Math(r"\frac{d}{dt}\mathbb{E}[" + step["source"] + r"] = \kappa\,\mathbb{E}[" + sp.latex(rhs) + "]"))

# %% [markdown]
# ## 5 · Initial states
#
# Start from one leaf and its exterior stem: **E=1, P1=P2=P3=0**. Integrate the generated linear expectation ODE.

# %%
KAPPA = 0.25
T = 3.

# %% [markdown]
# ## 6 · Dynamics

# %%
visual.close_controls(globals().get("controls"))
curves = visual.show_curve(visual.remy_curve(result, rate=KAPPA, end=T))
controls = visual.remy_sliders(result, curves, rate=KAPPA, end=T)

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
