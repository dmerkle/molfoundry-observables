# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Growing ring
#
# [← All examples](../paper_examples.ipynb)
#
# An unresolved closure search on rings of arbitrary size, with approximate dynamics. An additional example beyond the paper.

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
import inspection as files

# %% [markdown]
# ## 1 · Rules & state space
#
# Oriented rings of size ≥7, with three site states, conditional transitions, exchanges and state-copy growth.

# %%
from mfobs import analyze as analyze_model
from mfobs.applications import GrowingRing
backend = GrowingRing()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
rule_gmls = inputs.rule_gmls()
visual.rules([(n, rule_gmls[n]) for n in ('base_10', 'guard_10_not1', 'grow_copy_0')])
files.input_rules(rule_gmls)

# %% [markdown]
# ## 2 · Inputs & semantics
#
# One-site counts Q0,Q1,Q2. Their joint values do not determine all event rates.

# %%
visual.rules(list(inputs.query_gmls().items()), identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# Run 12 worklist expansions. Increase `MAX_EXPANSIONS` for a longer search.

# %%
MAX_EXPANSIONS = 12
result = analyze_model(inputs, worklist_limit=MAX_EXPANSIONS, include_evidence=True, backend=backend).record()
raw = result['linear_evidence']
display(view.overview(result))
suite.gml_table(raw['identity_patterns'])

# %% [markdown]
# ## 4 · Closure status
#
# Both closure searches remain unresolved. The dynamics below use a one-site approximation (R1) that factorizes correlations.

# %%
print("Linear closure:", result["linear"]["status"])
print("Polynomial jump closure:", result["pjc"]["status"])

# %% [markdown]
# ## 5 · Initial states
#
# Initial fractions for the approximate R1 model; these do not specify a unique graph state.

# %%
INITIAL = [.72, .14, .14]
visual.table([dict(zip(("q0","q1","q2"),INITIAL))])

# %% [markdown]
# ## 6 · Dynamics
#
# Compare g=0 with positive g. These curves show approximate site fractions.

# %%
T = 80.0
rates = {'r': 1.0, 'g': 0.12}
curves = suite.growing_curve(INITIAL, rates, end=T)
controls = suite.rate_controls(
    lambda rates, **opts: suite.growing_curve(INITIAL, rates, **opts), rates, curves['options'])

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
