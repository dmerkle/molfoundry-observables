"""Static graph, rule and query drawings through MolFoundry's public Python API."""
import molfoundry as mf
from mfobs.core.rules import rule_graph, is_identity
from mfobs.visual_links import compact_rule_gml, named_edge_labels


def graph_svg(graph):
    """Show the exact labelled graph, including explicit H and source vertex IDs."""
    printer = mf.GraphPrinter()
    printer.darkMode = mf.DarkMode.Off
    printer.drawMode = mf.DrawMode.Circles
    printer.withExternalId = True
    printer.withEdgeLabels = True
    printer.collapseHydrogens = False
    return graph.svg(printer)


def rule_svg(rule, *, identity=False):
    printer = mf.GraphPrinter()
    printer.darkMode = mf.DarkMode.Off
    printer.withIndex = True
    printer.withEdgeLabels = bool(named_edge_labels(rule.getGMLString()))
    printer.collapseHydrogens = False
    printer.drawRuleConditions = True
    if identity:
        if not is_identity(rule):
            raise ValueError("A single-panel drawing requires an identity rule")
        return rule_graph(rule).svg(printer)
    return rule.svg(printer)


def path_rule_svg(gml, *, identity=False):
    """Draw the incidence-encoded path with compact labels."""
    gml = compact_rule_gml(gml)
    rule = mf.Rule.fromGMLString(gml, add=False)
    if not identity:
        return rule_svg(rule)
    # P and N are graph types, not atoms: do not add implicit hydrogens.
    printer = mf.GraphPrinter()
    printer.drawMode = mf.DrawMode.Circles
    printer.darkMode = mf.DarkMode.Off
    printer.withIndex = True
    printer.withEdgeLabels = bool(named_edge_labels(gml))
    return rule_graph(rule).svg(printer)

