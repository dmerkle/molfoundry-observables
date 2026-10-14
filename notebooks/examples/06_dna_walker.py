# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # DNA walker
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/dna_walker_inspection.py) · [Model inputs](../../src/mfobs/models/dna_walker/inputs.py)
#
# Four graph rules yield equations for the extended and compressed walker configurations.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
NOTEBOOKS = PACKAGE / 'notebooks'
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
from IPython.display import display, Math
import visual_view as visual
import visual_suite as suite

# %% [markdown]
# ## 1 · Rules & state space
#
# One two-legged walker on a directed cycle of length ≥5. The legs attach to one site (compressed) or adjacent sites (extended). Auxiliary edge vertices with `src`/`tgt` links encode direction and distinguish the two legs.

# %%
from mfobs import analyze
from mfobs.applications import DNAWalker
backend = DNAWalker()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
visual.rules(list(inputs.rule_gmls().items()))

# %% [markdown]
# Model: Danos et al., *Rate Equations for Graphs*, CMSB 2020, §1.1.

# %% [markdown]
# ## 2 · Inputs & semantics
#
# The input query E counts extended configurations.

# %%
visual.rules(list(inputs.query_gmls().items()), identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# A short scaffold path extends uniquely from its leg attachments. The reducer uses this property to identify equivalent counts. The discovered coordinate O1 counts compressed configurations twice, once for each permutation of the two legs.

# %%
record = analyze(inputs, backend=backend).record()
print("Closure:", record["status"])
print("Discovered basis:", list(record["basis"]))
visual.rules([
    (name + " · " + record["display_names"][name],
     record["patterns"][next(iter(expression))])
    for name, expression in record["basis"].items() if name != "E"
], identity=True)

# %% [markdown]
# ## 4 · ODE
#
# Normalize the compressed count by two, then compute the stationary probabilities and mean velocity.

# %%
import sympy as sp
symbols = {n: sp.Symbol(n) for n in ("E", "C", "kFE", "kBE", "kFC", "kBC")}
comparison = record["comparison"]
for name, expression in comparison["physical_equations"].items():
    display(Math(r"\frac{d}{dt}\mathbb{E}[" + name + "] = " + sp.latex(sp.sympify(expression, locals=symbols))))
stationary = {symbols[n]: sp.sympify(value, locals=symbols)
              for n,value in comparison["stationary_probabilities"].items()}
velocity = sum(sp.Rational(row["displacement"]) * symbols[rate] * symbols[row["normalized_activity"]]
               for rate,row in comparison["movement_readout_from_rules"].items())
stationary_velocity = sp.factor(velocity.subs(stationary))
display(Math("v_{stationary} = " + sp.latex(stationary_velocity)))

# %% [markdown]
# ## 5 · Initial states
#
# One extended walker; probability E=1, C=0.

# %%
initial = {"E": 1., "C": 0.}
visual.table([initial])

# %% [markdown]
# ## 6 · Dynamics
#
# Integrate the computed two-dimensional equations. Counts are normalized to configuration probabilities: C = Craw/2.

# %%
T = 5.0
rates = {'kFE': 2.0, 'kBE': 0.4, 'kFC': 1.5, 'kBC': 0.3}
curves = suite.walker_curve(record, rates, end=T, initial=initial)
controls = suite.rate_controls(
    lambda rates, **opts: suite.walker_curve(record, rates, **opts), rates, {'end': T, 'initial': initial})

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
