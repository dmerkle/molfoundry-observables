"""Notebook presentation and numerical consumers for existing mfobs records.

No rule composition, observable discovery or model-specific moment equations
are implemented here. All mathematical inputs come from an analysis record.
Rule drawings use MolFoundry; application conditions are also listed separately.
"""
from __future__ import annotations

from html import escape
import re
from xml.dom.minidom import parseString

import molfoundry as mf
from mfobs.core.gml import condition_gml, constraint_gml
from mfobs.chemical_depiction import molecular_svg
from mfobs.visual_links import compact_rule_gml, editor_link, label_legend
from mfobs.rule_drawing import graph_svg, rule_svg

import pandas as pd
import sympy as sp
from IPython.display import HTML, Math, display
from mfobs.numerics.events import event_system, integrate_events, simulate_events
from mfobs.models.aba.analysis import word_coordinates


def analysis(record: dict) -> dict:
    return record.get("analysis", record)


def load_rule(gml: str):
    """Read a rule once through MolFoundry's public Python API."""
    return mf.Rule.fromGMLString(gml, add=False)


def rule_conditions_html(rule) -> str:
    """Present effective constraints and the public nested-condition AST."""
    labels = [f"NoEdge({c.source},{c.target})"
              for c in rule.matchConstraintsNoEdge if not c.implied]
    labels += [constraint_gml(c) for c in (
        *rule.matchConstraintsAdjacency, *rule.matchConstraintsLabelAny,
        *rule.matchConstraintsShortestPath)]
    condition = rule.matchCondition
    implied = [f"NoEdge({c.source},{c.target})"
               for c in rule.matchConstraintsNoEdge if c.implied]
    if implied:
        labels.append("created-edge admissibility (below)")
    if condition is not None:
        labels.append("nested application condition (below)")
    html = "<p><b>Application conditions:</b> " + (
        " · ".join(f"<code>{escape(label)}</code>" for label in labels) or "none"
    ) + "</p>"
    if implied:
        html += ("<details><summary>Created-edge admissibility</summary><p>"
                 + escape(" · ".join(implied))
                 + "</p><p>These edges must be absent before creation. "
                   "The conditions may also be written explicitly in the input GML.</p></details>")
    if condition is not None:
        html += ("<details><summary>Nested application condition</summary><pre>"
                 + escape(condition_gml(condition)) + "</pre></details>")
    return html


def notebook_svg(drawing: str) -> str:
    """Keep MolFoundry's SVG colours when notebook previews discard stylesheets."""
    root = parseString(drawing).documentElement
    _, _, width, height = root.getAttribute("viewBox").split()
    root.setAttribute("width", width)
    root.setAttribute("height", height)
    for style in list(root.getElementsByTagName("style")):
        # MolFoundry uses tag selectors and SVG-id-scoped class selectors.
        elements = list(style.parentNode.getElementsByTagName("*"))
        for selector, declarations in re.findall(r"([^{}]+)\{([^{}]+)\}",
                                                  style.firstChild.data):
            selector = selector.strip().split()[-1]
            attributes = [item.split(":", 1) for item in declarations.split(";")
                          if item.strip()]
            for element in elements:
                matches = (selector[1:] in element.getAttribute("class").split()
                           if selector.startswith(".") else element.tagName == selector)
                if matches:
                    for key, value in attributes:
                        element.setAttribute(key.strip(), value.strip())
        style.parentNode.removeChild(style)
    # GitHub also turns literal UTF-8 characters into visible byte escapes.
    return root.toxml(encoding="ascii").decode("ascii")


def show_molfoundry_rule(gml: str, title: str = "DPO rule", *, identity=False, legend=True) -> None:
    """Display the public MolFoundry SVG, conditions and an editor link."""
    shown_gml = compact_rule_gml(gml)
    rule = load_rule(shown_gml)
    display(HTML(f"<h4>{escape(title)}</h4>"))
    if shown_gml != gml:
        from mfobs.rule_drawing import path_rule_svg
        drawing = path_rule_svg(gml, identity=identity)
    else:
        drawing = rule_svg(rule, identity=identity)
    display({"image/svg+xml": notebook_svg(drawing)}, raw=True)
    if legend and label_legend(gml):
        display(HTML("<p>" + escape(label_legend(gml)) + "</p>"))
    display(HTML(rule_conditions_html(rule) + "<p>" + editor_link(gml) + "</p>"))


def show_graph(gml: str, title="Host graph") -> None:
    graph = mf.Graph.fromGMLString(gml, add=False)
    display(HTML(f"<h4>{escape(title)}</h4>"))
    display({"image/svg+xml": notebook_svg(graph_svg(graph))}, raw=True)


def show_molecule(gml: str, title="Molecule", *, indices: bool = False) -> None:
    """Draw the native graph with MolFoundry, retaining source IDs when requested."""
    graph = mf.Graph.fromGMLString(gml, add=False)
    display({"image/svg+xml": notebook_svg(molecular_svg(
        graph, title, width=760, height=240, indices=indices))}, raw=True)


def show_events(record: dict) -> pd.DataFrame:
    events = event_system(record)
    rows = []
    for name, rate, delta, activity in zip(events.identifiers, events.channels, events.deltas, events.activities):
        rows.append({"Event": name, **{f"Δ {n}": int(v) for n, v in zip(events.names, delta)}, "Rate": f"{rate} × ({activity})"})
    return pd.DataFrame(rows, columns=["Event", *(f"Δ {n}" for n in events.names), "Rate"]).set_index("Event")


def show_equations(record: dict, names=None) -> None:
    events = event_system(record)
    result = analysis(record)
    closed = result.get("linear", {}).get("status") == "CLOSED"
    display(HTML("<p><b>" + ("Closed expectation equations." if closed else "First-moment identities; linear closure unresolved.") + "</b></p>"))
    if names is not None and set(names) - set(events.names):
        raise ValueError("Unknown requested equation coordinate")
    for name, expr in zip(events.names, events.drift()):
        if names is not None and name not in names:
            continue
        poly = sp.Poly(expr, *events.symbols)
        parts = []
        for exponents, coefficient in poly.terms():
            if not coefficient:
                continue
            monomial = sp.prod(s**p for s, p in zip(events.symbols, exponents))
            expectation = sp.latex(monomial) if monomial == 1 else r"\mathbb E\!\left[" + sp.latex(monomial) + r"\right]"
            parts.append("(" + sp.latex(coefficient) + ")" + expectation)
        display(Math(r"\frac{d}{dt}\mathbb E[" + sp.latex(sp.Symbol(name)) + "]=" + (" + ".join(parts) or "0")))


def overview(record: dict) -> pd.DataFrame:
    result = analysis(record)
    return pd.DataFrame([{
        "Linear closure": result.get("linear", {}).get("status", "not run"),
        "Event model": result.get("pjc", {}).get("status", "not run"),
        "Coordinates": len(result["pjc"]["observable_coordinates"]) if result.get("pjc", {}).get("observable_coordinates") is not None else None,
        "Event rows": len(result["pjc"]["operator"].get("terms", [])) if result.get("pjc", {}).get("operator") is not None else None,
    }])
