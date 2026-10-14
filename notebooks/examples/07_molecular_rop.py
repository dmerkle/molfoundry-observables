# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Molecular ring-opening polymerization
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/rop_analysis.py) · [Model inputs](../../src/mfobs/models/molecular_rop/inputs.py)
#
# Discover counts and events from atom-resolved lactone rules. The default query counts ABA in an unbranched polyester chain.

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
# A = δ-valerolactone; B = ε-caprolactone. One benzyl-initiated chain, free lactones and an optional gate. The graph family excludes branching, backbiting and transesterification.

# %%
from mfobs.applications import MolecularROP
from mfobs.models.molecular_rop import inputs as rop
from mfobs.core.gml import graph_gml
CASE = 'aba'
MOTIF = 'ABA'  # Also try 'ABBA'.
if CASE in ('baseline', 'noedge'):
    backend = MolecularROP(CASE)
elif CASE in ('aba', 'aba-gated'):
    backend = MolecularROP.motif(MOTIF, gated=CASE.endswith('gated'))
else:
    raise ValueError('Unknown molecular case')
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
from mfobs import ModelInputs
inputs = ModelInputs(Omega, rules, observables)
selection = backend.configuration()
rule_gmls = inputs.rule_gmls()
visual.rules([('insert_' + selection.alphabet[0], rule_gmls['insert_' + selection.alphabet[0]])])
files.input_rules(rule_gmls)
for symbol in selection.alphabet:
    view.show_molecule(graph_gml(rop.monomer(symbol)), 'Free lactone ' + symbol)

# %% [markdown]
# ## 2 · Inputs & semantics
#
# The ABA cases query the ABA motif. The B-only cases query the number of B units.

# %%
visual.rules(list(inputs.query_gmls().items()), identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# Compose the molecular rules and reduce counts on the specified chain family.

# %% [markdown]
# **Algorithm 1 · Linear expectation closure**

# %%
linear = backend.linear_closure(inputs, maximum_dimension=24)
print("Linear expectation closure:", linear["status"])

# %% [markdown]
# **Algorithm 2 · Pattern support and polynomial events**

# %%
result = backend.jump_closure(inputs)
visual.support(result)
suite.gml_table(result['support_discovery']['patterns'])

# %% [markdown]
# ## 4 · pJC
#
# The event model closes. A finite event model does not by itself give closed equations for expectations.

# %%
result["linear"] = linear
display(view.overview(result))
display(view.show_events(result))
view.show_equations(result, names=result['support_discovery']['selected_queries'])

# %% [markdown]
# ## 5 · Initial states
#
# Choose a chain and the free monomer pools.

# %%
WORD = "BBBB" if CASE in ("baseline", "noedge") else "ABAB"
POOLS = {s: (2 if s == "A" else 1) for s in selection.alphabet}
OPEN = True
view.show_molecule(graph_gml(rop.chain(tuple(WORD))), "Benzyl-initiated chain " + WORD)
initial = view.word_coordinates(result, WORD, POOLS, opened=OPEN)
visual.table([initial])

# %% [markdown]
# ## 6 · Dynamics
#
# Mean and variance from finite-state projection (FSP), compared with SSA estimates. FSP truncation bounds and sampling uncertainty are reported below.

# %%
T = 4.0
RUNS = 100
SEED = 20260924
rates = {'beta_A': 2.0, 'beta_B': 1.5, 'k_A': 0.9, 'k_B': 0.01}
rates = {k: v for k, v in rates.items() if k.rsplit('_', 1)[-1] in selection.alphabet}
if selection.gated:
    rates.update(a=0.4, b=0.6)
curves = visual.show_curve(visual.aba_curve(result, initial, rates, end=T, paths=RUNS, seed=SEED))
controls = visual.aba_sliders(result, initial, rates, curves)

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
