# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Composition measurements
#
# [← All examples](../paper_examples.ipynb)
#
# Inspect recorded MolFoundry and FRAGGER timings for the same composition workload. This notebook does not rerun the external benchmark.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
NOTEBOOKS = PACKAGE / 'notebooks'
sys.path[:0] = [str(NOTEBOOKS), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import pandas as pd
from IPython.display import display

# %% [markdown]
# ## 1 · Rules & state space
#
# Six terminal-extension rules: p→T becomes p→s→T, with p∈{R,A,B}, s∈{A,B}. No vertices are deleted. FRAGGER uses directed multigraphs; MolFoundry uses incidence graphs.

# %%
EVIDENCE = PACKAGE / 'notebooks/visual/evidence/ab-scaling'

# %% [markdown]
# ## 2 · Inputs & semantics
#
# The query is an alternating AB word of length n. Each measurement composes it in both orders with all six rules: twelve calls.

# %%
summary = pd.read_csv(EVIDENCE / "summary.csv")

# %% [markdown]
# ## Composition measurements
#
# MolFoundry 0.1.9 (14 October) and FRAGGER (9 October), both on `algochem-compute` (Linux/Xeon w7-2595X). FRAGGER samples are retained unchanged. Medians of five runs after two per-case warm-ups; startup and compilation excluded. [Protocol](https://stochmech.algochem.techfak.de/experiments/comparison/ab_scaling/README.md).

# %%
timings = summary[["n", "molfoundry_median_seconds", "fragger_median_seconds"]].rename(
    columns={"molfoundry_median_seconds": "MolFoundry (s)",
             "fragger_median_seconds": "FRAGGER (s)"})
display(timings.style.hide(axis="index").format(
    {"MolFoundry (s)": "{:.3g}", "FRAGGER (s)": "{:.3g}"}))

# %% [markdown]
# FRAGGER is faster at some small sizes. These are composition timings, not complete closure runs.

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
