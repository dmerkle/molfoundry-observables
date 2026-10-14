# Generated from the adjacent .ipynb by scripts/notebook_sources.py.
# Edit the notebook, then regenerate this file.

# %% [markdown]
# # Chlorination with a neighbour NAC
#
# [← All examples](../paper_examples.ipynb) · [Python analysis](../../tutorials/chlorination_analysis.py) · [Model inputs](../../src/mfobs/models/chlorination/inputs.py)
#
# Compare the same chlorination rule with and without its neighbour condition. The graph family has an arbitrary finite alkyl tail; the generated event model closes in both variants.

# %%
from pathlib import Path
import sys
PACKAGE = next((p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'src/mfobs').is_dir()))
sys.path[:0] = [str(PACKAGE / 'notebooks'), str(PACKAGE / 'src'), str(PACKAGE)]
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import visual_view as visual

# %% [markdown]
# ## Rules and query
#
# All five ring bonds stay aromatic (`:`). The two chlorination rules differ only in the NAC. The query Z counts ring chlorines.

# %%
from mfobs import ModelInputs, analyze
from mfobs.applications import Chlorination

backends = {v: Chlorination(v) for v in ("guarded", "ablated")}
defaults = {v: backend.specification() for v, backend in backends.items()}
Omega = defaults["guarded"].omega
rules = {v: model.rules for v, model in defaults.items()}
observables = defaults["guarded"].observables
inputs = {v: ModelInputs(Omega, rules[v], observables) for v in backends}
from mfobs.visual_links import editor_url
from IPython.display import Markdown, display

links = []
for variant, model in inputs.items():
    for name, gml in model.rule_gmls().items():
        if variant == "ablated" and name == "grow":
            continue
        label = "Tail growth" if name == "grow" else (
            "With NAC" if variant == "guarded" else "Without NAC")
        links.append(f"[{label}]({editor_url(gml)})")
for name, gml in inputs["guarded"].query_gmls().items():
    links.append(f"[Query {name}]({editor_url(gml)})")
display(Markdown(" · ".join(links)))

# %% [markdown]
# ## 3 · Observable discovery

# %%
results = {v: analyze(model, backend=backends[v]).record() for v, model in inputs.items()}
visual.table([
    {"Rule": "With NAC" if v == "guarded" else "Without NAC",
     "Counts": ", ".join(r["pjc"]["coordinates"]),
     "pJC": r["pjc"]["status"], "Events": len(r["events"])}
    for v, r in results.items()
])

# %% [markdown]
# ## 4 · pJC

# %%
for variant, result in results.items():
    print("With NAC" if variant == "guarded" else "Without NAC")
    visual.table([{"Rate / κ": event["activity"]["text"],
                   "Jump": ", ".join(f"Δ{name}={value:+d}" for name, value in event["delta"].items())}
                  for event in result["events"]])

# %% [markdown]
# ## 5 · Initial states
#
# `INITIAL_MASK=0` starts unchlorinated; `1` starts with one chlorine, `9` with chlorine at sites 1 and 4. The curves integrate the finite forward equations of the generated event model. Tail growth has zero jump on these counts.

# %%
INITIAL_MASK = 0  # also try 1 or 9
RATE = 1.
T = 3.

# %% [markdown]
# ## 6 · Dynamics

# %%
visual.close_controls(globals().get("controls"))
curves = visual.show_curve(visual.context_nac_curve(results, mask=INITIAL_MASK, rate=RATE, end=T))
controls = visual.context_nac_sliders(results, curves, mask=INITIAL_MASK, rate=RATE, end=T)

# %% [markdown]
# [← All examples](../paper_examples.ipynb)
