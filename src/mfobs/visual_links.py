"""Compact labels and MolFoundry inspection links for notebook drawings."""
from html import escape
import re
from urllib.parse import quote

# Keep downloads for graph URLs exceeding the current HTTP/2 endpoint limit.
MAX_GRAPH_LINK_LENGTH = 4000

# Only display copies are renamed; analysis inputs retain their original labels.
DISPLAY_LABELS = {
    "Position": ("P", "position"), "Next": ("N", "next-link node"),
    "Terminal": ("T", "terminal-link node"), "End": ("E", "end"),
    "Root": ("R", "root"), "Symbol:A": ("A", "letter A"),
    "Symbol:B": ("B", "letter B"), "Feed:A": ("FA", "A feed"),
    "Feed:B": ("FB", "B feed"), "Feed-A": ("FA", "A feed"),
    "Feed-B": ("FB", "B feed"), "GateSource": ("GS", "gate source"),
    "GateTarget": ("GT", "gate target"), "bimotor": ("M", "motor"),
    "chain": ("C", "scaffold site"), "edge_dna": ("D", "scaffold-link node"),
    "edge_leg": ("L", "leg-link node"),
}


def compact_rule_gml(gml):
    """Abbreviate graph-type labels, also inside application conditions."""
    return re.sub(r'\blabel\s+"([^"\\]*(?:\\.[^"\\]*)*)"',
                  lambda m: 'label "' + DISPLAY_LABELS[m[1]][0] + '"'
                  if m[1] in DISPLAY_LABELS else m[0], gml)


def label_legend(gml):
    """One short legend for the graph types occurring in these drawings."""
    labels = set(re.findall(r'\blabel\s+"([^"\\]*)"', gml))
    meanings = dict(DISPLAY_LABELS[label] for label in DISPLAY_LABELS if label in labels)
    return "; ".join(f"{short} = {meaning}" for short, meaning in meanings.items())


def named_edge_labels(gml):
    """Edge types that need text labels rather than chemical bond strokes."""
    labels = re.findall(r'\bedge\s*\[[^\[\]]*\blabel\s+"([^"\\]*)"', gml)
    return set(labels) - {"", "-", "=", "#", ":"}


def editor_url(gml):
    """Open the abbreviated drawing, showing structural edge types."""
    shown = compact_rule_gml(gml)
    if named_edge_labels(gml):
        url = ("https://molfoundry.de/depiction?gml=" + quote(shown, safe="")
               + "&drawMode=circles&withEdgeLabels=on")
        if re.match(r"\s*rule\s*\[", shown):
            url += "&drawRuleConditions=on"
        return url
    if re.match(r"\s*graph\s*\[", shown):
        return "https://molfoundry.de/depiction?gml=" + quote(shown, safe="") + "&drawMode=circles"
    return "https://molfoundry.de/rule-editor?gml=" + quote(shown, safe="")


def editor_link(gml, label="View in MolFoundry"):
    return (f'<a href="{escape(editor_url(gml), quote=True)}" '
            f'target="_blank" rel="noopener">{escape(label)}</a>')

