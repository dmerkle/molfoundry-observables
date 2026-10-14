# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Fixed ring
#
# [← All examples](../../paper_examples.ipynb)
#
# The eight-site binding model restricted to a single ring containing all eight sites. An additional example beyond the paper.
#
# [State space and dynamics](https://stochmech.algochem.techfak.de/experiments/ring8/index.html)

# %%
from pathlib import Path
import sys
PACKAGE = next(p for p in (Path.cwd(), *Path.cwd().parents)
               if (p / "src/mfobs").is_dir())
sys.path.insert(0, str(PACKAGE / "src"))
from mfobs import Ring8, analyze
backend = Ring8()
inputs = backend.specification()
Omega = inputs.omega
rules = inputs.rules
observables = inputs.observables

# %% [markdown]
# ## 2 · Inputs & semantics
#
# Rules and queries are read from the same model manifest used for analysis.

# %%
from mfobs.visual_links import editor_url
from IPython.display import Markdown, display

for heading, gmls in (("Rules", inputs.rule_gmls()), ("Queries", inputs.query_gmls())):
    links = [f"[{name}]({editor_url(gml)})" for name, gml in gmls.items()]
    display(Markdown(heading + ": " + " · ".join(links)))
print("Selected queries:", list(observables))

# %% [markdown]
# ## 3 · Observable discovery

# %%
result = analyze(inputs, include_evidence=True, backend=backend).record()
linear = result["linear_evidence"]
print("Linear closure:", linear["status"])
print("Discovered basis:", linear["basis"])
print("Polynomial jump closure:", result["pjc"]["status"])
print("Event classes:", result["pjc"]["event_class_count"])

# %% [markdown]
# ## 4 · ODE
#
# The equations below use the computed coefficients. The complete generated event model is in `result["pjc"]`.

# %%
from html import escape
from IPython.display import HTML
from mfobs.outputs import ode_text

equations = [line for line in ode_text(linear).splitlines() if line.startswith("dmu_")]
display(HTML("<pre>" + escape(equations[0]) + "</pre>"
             "<details><summary>Remaining equations</summary><pre>"
             + escape("\n".join(equations[1:])) + "</pre></details>"))

# %% [markdown]
# [← All examples](../../paper_examples.ipynb)
