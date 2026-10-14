# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # ABA · change the query
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/query_variation_analysis.py) · [Model inputs](../../src/mfobs/models/aba/inputs.py)
#
# Compare discovery for five motifs with the same rules and search budgets.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
NOTEBOOKS = PACKAGE / 'notebooks'
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import paper_view as view
import visual_view as visual
import visual_suite as suite
import inspection as files

# %% [markdown]
# ## 1 · Rules & state space
#
# Four ungated feed/append rules are shared by all five queries.
#
# Directed links use `src`/`tgt` edges.

# %%
from mfobs import ModelInputs, ABA
MOTIFS = ("AB", "ABA", "ABBA", "AAAA", "ABABA")
backends = {motif: ABA.motif(motif, gated=False) for motif in MOTIFS}
defaults = {motif: backend.specification() for motif, backend in backends.items()}
Omega = defaults[MOTIFS[0]].omega
rules = defaults[MOTIFS[0]].rules
observables = {motif: model.observables for motif, model in defaults.items()}
inputs = {motif: ModelInputs(Omega, rules, observables[motif]) for motif in MOTIFS}
rule_gmls = inputs[MOTIFS[0]].rule_gmls()
visual.rules([("insert_A", rule_gmls["insert_A"])])
files.input_rules(rule_gmls)

# %% [markdown]
# ## 2 · Inputs & semantics
#
# Vary only the motif. Pattern and expansion budgets are 64; the linear dimension budget is 24.

# %%
visual.rules([(name, gml) for model in inputs.values()
              for name, gml in model.query_gmls().items()], identity=True)

# %% [markdown]
# ## 3 · Observable discovery
#
# Run the two searches independently from the same rules and query.

# %% [markdown]
# **Algorithm 1 · Linear expectation closure**

# %%
linear = {motif: backends[motif].linear_closure(model, maximum_dimension=24)
          for motif, model in inputs.items()}

# %% [markdown]
# **Algorithm 2 · Pattern support and polynomial events**

# %%
results = {motif: backends[motif].jump_closure(
    model, support_pattern_limit=64, support_expansion_limit=64)
    for motif, model in inputs.items()}

# %%
for motif, result in results.items():
    result["linear"] = linear[motif]
visual.table([
    {"Query": motif, "Counts": ", ".join(r["support_discovery"]["patterns"]),
     "pJC": r["pjc"]["status"], "Linear closure": r["linear"]["status"]}
    for motif, r in results.items()
])
for motif, result in results.items():
    suite.gml_table(result["support_discovery"]["patterns"], motif + " · patterns")

# %% [markdown]
# ## 4 · Generator identities
#
# Each row uses its own generated query drift and event model.

# %%
for motif,r in results.items():
    print(motif)
    view.show_equations(r, names=r["support_discovery"]["selected_queries"])

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
