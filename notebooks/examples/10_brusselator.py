# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Brusselator–ROP
#
# [← All examples](../paper_examples.ipynb)
#
# Exact events and stochastic simulation, compared with approximate oscillatory dynamics. An additional example beyond the paper.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
NOTEBOOKS = PACKAGE / 'notebooks'
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
from IPython.display import display
import paper_view as view
import visual_view as visual
import visual_suite as suite

# %% [markdown]
# ## 1 · Rules & state space
#
# Six rules implement the Brusselator reactions, monomer input and X-catalysed lactone growth.

# %%
from mfobs.applications import BrusselatorROP
backend = BrusselatorROP()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
visual.rules(list(inputs.rule_gmls().items()))

# %% [markdown]
# ## 2 · Inputs & semantics
#
# Start with population counts X and Y; discovery adds monomer and chain-length counts.
#
#
# Molecule counts divide pattern matches by molecular automorphisms.

# %%
visual.rules(list(inputs.query_gmls().items()), identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# Run linear closure and event construction independently from the rules.

# %% [markdown]
# **Algorithm 1 · Linear expectation closure**

# %%
linear = backend.linear_closure(inputs, maximum_dimension=24)
print("Linear expectation closure:", linear["status"])

# %% [markdown]
# **Polynomial events**

# %%
result = backend.jump_closure(inputs)

# %%
result["linear"] = linear
display(view.overview(result))

# %% [markdown]
# ## 4 · pJC
#
# The event model is exact. First-moment equations require higher moments; the deterministic approximation factorizes them.

# %%
view.show_equations(result)
display(view.show_events(result))

# %% [markdown]
# ## 5 · Initial states
#
# Initial concentrations x,y and counts M,N. X=Vx and Y=Vy for SSA; volume V is in the dynamics settings.

# %%
INITIAL = [.4, 2.5, 4., 0.]
visual.table([dict(zip(("x","y","M","N"),INITIAL))])

# %% [markdown]
# ## 6 · Dynamics
#
# Compare b=1.8 (damped) and b=2.2 (oscillatory) using the rate control. SSA uses the computed falling-factorial event activity. The band is uncertainty of the sampled mean, not a certified expectation curve.

# %%
T = 60.0
RUNS = 32
SEED = 20260924
VOLUME = 40
rates = dict(a=1.0, b=2.2, k1=1.0, k2=1.0, k3=1.0, k4=1.0, beta=0.3, k_g=0.15)
curves = suite.bruss_curve(result, INITIAL, rates, end=T, paths=RUNS, seed=SEED, volume=VOLUME)
controls = suite.rate_controls(
    lambda rates, **opts: suite.bruss_curve(result, INITIAL, rates, **opts), rates, curves['options'])

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
