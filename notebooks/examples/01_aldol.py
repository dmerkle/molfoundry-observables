# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Aldol addition
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/aldol_analysis.py) · [Model inputs](../../examples/tcs/aldol.gml)
#
# One NoEdge condition changes which compositions exist and which signed terms cancel.

# %%
from pathlib import Path
import sys
NOTEBOOKS = next(p for p in (Path.cwd(), *Path.cwd().parents)
                 if (p / "paper_view.py").is_file())
PACKAGE = NOTEBOOKS.parent
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / "src"), str(PACKAGE)]
import matplotlib
matplotlib.use("module://matplotlib_inline.backend_inline")
import visual_view as visual

EXAMPLE = 'visual-aldol'
run = visual.Run(EXAMPLE)
import inspection as files

# %% [markdown]
# ## Model
#
# The aldol rule retains its NoEdge condition. We change only the query condition. A valence filter checks atom labels and bond-order sums (C ≤ 4, O ≤ 2, H ≤ 1). This example studies local compositions, rather than a full chemical state space.

# %%
from mfobs.state_space import covalent_compositions
from examples.tcs.replay import load_rules
Omega = covalent_compositions()
inputs = load_rules(omega=Omega)
visual.rules([('Aldol addition · L ← K → R', inputs['aldol'])])

# %% [markdown]
# Blue: deleted · Green: created · Red ×: NoEdge. Vertex numbers identify the atoms in the GML.

# %% [markdown]
# ## Query
#
# Count the same pattern with and without **NoEdge(0, 5)**. Each picture shows the input of its identity rule.

# %%
visual.rules([
    ("Without query NoEdge", inputs["query_unconditioned"]),
    ("With query NoEdge", inputs["query_conditioned"]),
], identity=True)

# %% [markdown]
# ## Composition
#
# Compose in both orders: reaction → query counts **after**; query → reaction counts **before**. The signed contribution is after minus before.

# %%
from examples.tcs.compare_conditions import compare_conditions
with run.capture('compositions'):
    comparison = compare_conditions(run.path / 'compositions', export_candidates=True)
from examples.tcs.table1 import export_computed_table
with run.capture('plain-table'):
    export_computed_table(run.path / 'plain-table')

# %%
visual.table([
    {"Query NoEdge": label,
     "After compositions": comparison["variants"][key]["orders"]["reaction_then_query"]["filtered"],
     "Before compositions": comparison["variants"][key]["orders"]["query_then_reaction"]["filtered"],
     "Nonzero signed terms": comparison["variants"][key]["nonzero_signed_terms"]}
    for label, key in [("absent", "reaction_noedge"), ("present", "reaction_and_query_noedge")]
], index="Query NoEdge")

# %%
aldol_files = files.aldol_tables(
    comparison, run.path / "compositions", run.path / "plain-table", run)

# %% [markdown]
# The query condition excludes three after-count classes (`r68`, `r73`, `r120`). It also prevents cancellations. Below is the selected `r180` pair: its positive graph spans coincide, but its application conditions differ.

# %%
pair = comparison["comparison"]["selected_r180"]["reaction_and_query_noedge"]
visual.rules([(("after" if row["after"] else "before") + "-count composite",
              (run.path / "compositions" / row["gml"]).read_text(encoding='utf-8'))
              for row in pair])

# %% [markdown]
# Source: Behr et al., *Rewriting Theory for the Life Sciences*, §8.3, [doi:10.1016/j.tcs.2021.07.026](https://doi.org/10.1016/j.tcs.2021.07.026).

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
