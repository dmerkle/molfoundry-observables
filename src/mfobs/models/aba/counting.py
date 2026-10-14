"""Exact counts of composed path patterns on Ω: pools, motifs, suffixes and the gate."""
from __future__ import annotations
from collections import Counter
import networkx as nx
import sympy as sp
import molfoundry as mf
from mfobs.core.engine import commutator
from mfobs.core.pattern_support import Image, Pattern, UnsupportedReduction, discover_support, reduce_signed_rules
from mfobs.models.aba.inputs import SourceRule, graph_rule, query_pattern, source_rules, state_family


def query_motifs(selection) -> tuple[str, ...]:
    """Recognise common file-based queries as canonical path motifs."""

    motifs = []
    for query in selection.observables:
        if len(query.terms) != 1 or query.terms[0].coefficient != 1:
            raise ValueError("ROP accepts one unit-coefficient motif identity per query")
        rule = mf.Rule.fromGMLString(query.terms[0].gml, add=False)
        image = reduce_input(rule)
        if len(image.factors) != 1:
            raise ValueError(f"{query.name}: GML is not one terminal-path motif")
        pattern = next(iter(image.factors.values()))
        if pattern.kind != "motif" or image.expression != sp.Symbol(pattern.name):
            raise ValueError(f"{query.name}: GML is not one canonical motif identity")
        if query.name != pattern.name:
            raise ValueError(
                f"{query.name}: canonical name for this GML is {pattern.name}"
            )
        motifs.append(pattern.metadata["word"])
    return tuple(motifs)


def side_graph(rule, side="left") -> nx.Graph:
    """Read a native rule side directly for directed-path recognition."""
    if side not in ("left", "right"):
        raise ValueError(f"Unknown rule side {side!r}")
    label = side + "Label"
    graph = nx.Graph()
    graph.add_nodes_from((vertex.id, {"label": value}) for vertex in rule.vertices
                         if (value := getattr(vertex, label)) is not None)
    graph.add_edges_from((edge.source.id, edge.target.id, {"label": value})
                         for edge in rule.edges
                         if (value := getattr(edge, label)) is not None)
    return graph


def structural_violation(rule, side="left") -> str | None:
    graph = side_graph(rule, side)
    labels = nx.get_node_attributes(graph, "label")
    counts = Counter(labels.values())
    for label, number in counts.items():
        if (label in ("End", "Root", "GateSource", "GateTarget") or label.startswith("Feed:")) and number > 1:
            return "unique-language-object:" + label
    directed = nx.DiGraph()
    for i, label in labels.items():
        neighbours = list(graph[i])
        if label in ("Next", "Terminal"):
            src = [j for j in neighbours if graph[i][j]["label"] == "src"]
            tgt = [j for j in neighbours if graph[i][j]["label"] == "tgt"]
            if len(src) > 1 or len(tgt) > 1 or len(src) + len(tgt) != len(neighbours):
                return "nonfunctional-incidence"
            if src and labels[src[0]] != "Position":
                return "source-type"
            if tgt and labels[tgt[0]] != ("Position" if label == "Next" else "End"):
                return "target-type"
            if src and tgt:
                directed.add_edge(src[0], tgt[0])
        elif label == "Position":
            for bond in ("kind", "src", "tgt"):
                if sum(graph[i][j]["label"] == bond for j in neighbours) > 1:
                    return "position-functionality:" + bond
        elif label.startswith("Symbol:") or label in ("Root", "End"):
            if len(neighbours) > 1:
                return "unique-owner"
    if not nx.is_directed_acyclic_graph(directed):
        return "directed-cycle"
    return None


def projected_input(rule):
    """Retain conditions and check when the full-input diagonal is admissible."""
    graph = side_graph(rule)
    known = {"Position", "Root", "End", "Next", "Terminal", "GateSource", "GateTarget"}
    for _, label in graph.nodes(data="label"):
        if label not in known and not label.startswith(("Symbol:", "Feed:")):
            raise UnsupportedReduction("node outside the supported terminal-path signature")
    for left, right, bond in graph.edges(data="label"):
        labels = {graph.nodes[left]["label"], graph.nodes[right]["label"]}
        valid = (
            (bond == "kind" and "Position" in labels and
             any(label == "Root" or label.startswith("Symbol:") for label in labels))
            or (bond == "src" and labels in ({"Position", "Next"}, {"Position", "Terminal"}))
            or (bond == "tgt" and labels in ({"Position", "Next"}, {"End", "Terminal"}))
            or (bond == "bound" and labels == {"GateSource", "GateTarget"})
        )
        if not valid:
            raise UnsupportedReduction("bond outside the supported terminal-path signature")
    if rule.matchCondition is not None:
        raise UnsupportedReduction("generic nested matchCondition needs a supported reducer")
    if rule.matchConstraintsLabelAny or rule.matchConstraintsShortestPath:
        raise UnsupportedReduction("unhandled label/path condition")
    for vertex in rule.vertices:
        if vertex.membership.value != "deleted":
            continue
        i = vertex.id
        if (graph.nodes[i]["label"] not in ("Next", "Terminal") or len(graph[i]) != 2
                or {graph[i][j]["label"] for j in graph[i]} != {"src", "tgt"}):
            raise UnsupportedReduction("DPO interface cannot be erased: unsaturated deleted vertex")
    no_edges = sorted(set(tuple(sorted((c.source, c.target)))
                          for c in rule.matchConstraintsNoEdge))
    if any(graph.has_edge(i, j) for i, j in no_edges):
        return None
    adjacency = []
    for condition in rule.matchConstraintsAdjacency:
        i, op, value = condition.id, condition.op.value, condition.count
        if graph.nodes[i]["label"] != "Position":
            raise UnsupportedReduction("adjacency condition outside supported Position counts")
        ns, es = condition.nodeLabels or (), condition.edgeLabels or ()
        if set(ns) != {"Next", "Terminal"} or set(es) != {"src", "tgt"} or op not in ("=", ">") or value != 0:
            raise UnsupportedReduction("unhandled adjacency condition: " + str(condition))
        number = sum(graph.nodes[j]["label"] in ns and graph[i][j]["label"] in es for j in graph[i])
        if op == "=" and number > 0:
            return None
        if op == ">" and number > 0:
            continue
        adjacency.append((i, op, value))
    adjacency = sorted(set(adjacency))
    for i, op, value in adjacency:
        if (i, ">" if op == "=" else "=", value) in adjacency:
            return None
    return graph, adjacency, no_edges


def recognize_component(graph, adjacency, no_edges):
    """Name a graph that was already produced; never enumerate candidate motifs."""
    labels = nx.get_node_attributes(graph, "label")
    edges = [(i, j, d["label"]) for i, j, d in graph.edges(data=True)]
    identity = graph_rule(labels, edges, adjacency=adjacency, no_edges=no_edges)
    if len(graph) == 1 and (next(iter(labels.values())).startswith("Feed:") or
                            next(iter(labels.values())) in ("GateSource", "GateTarget", "End")):
        return sp.Integer(1), None, "unique-object"
    if len(graph) == 2 and set(labels.values()) == {"GateSource", "GateTarget"}:
        if edges and (len(edges) != 1 or edges[0][2] != "bound"):
            raise UnsupportedReduction("unexpected gate bond")
        if not edges and not no_edges:
            return sp.Integer(1), None, "unique-gate-pair"
        opened = graph_rule(labels, [], no_edges=[tuple(labels)])
        pattern = Pattern("X", opened, "gate", {})
        return (1 - sp.Symbol("X") if edges else sp.Symbol("X")), pattern, "unique-gate"
    if no_edges:
        raise UnsupportedReduction("cross-component absence condition is not separable here")
    positions = [i for i, label in labels.items() if label == "Position"]
    directed = nx.DiGraph()
    directed.add_nodes_from(positions)
    tags, terminal = {}, None
    for i, label in labels.items():
        if label == "Position":
            attached = [labels[j][7:] for j in graph[i]
                        if graph[i][j]["label"] == "kind" and labels[j].startswith("Symbol:")]
            if attached:
                tags[i] = attached[0]
        elif label in ("Next", "Terminal"):
            src = [j for j in graph[i] if graph[i][j]["label"] == "src"]
            tgt = [j for j in graph[i] if graph[i][j]["label"] == "tgt"]
            if len(src) != 1 or len(tgt) != 1:
                raise UnsupportedReduction("partial incidence edge")
            if label == "Next":
                directed.add_edge(src[0], tgt[0])
            else:
                terminal = src[0]
    if len(positions) == 1 and len(graph) == 2 and len(edges) == 1 and edges[0][2] == "kind":
        symbol = tags.get(positions[0])
        if symbol is None:
            raise UnsupportedReduction("untyped token count")
        if adjacency == [(positions[0], "=", 0)]:
            name, kind = "M_" + symbol, "pool"
        elif adjacency == [(positions[0], ">", 0)]:
            name, kind = "N_" + symbol, "motif"
        elif not adjacency:
            name, kind = "T_" + symbol, "total"
        else:
            raise UnsupportedReduction("unsupported token condition")
        return sp.Symbol(name), Pattern(name, identity, kind, {"word": symbol}), kind
    if not positions or not nx.is_weakly_connected(directed):
        raise UnsupportedReduction("component outside the connected path-count language")
    ordered = list(nx.topological_sort(directed))
    word = "".join(tags.get(i, "?") for i in ordered)
    if terminal is not None and ordered[-1] != terminal:
        raise UnsupportedReduction("terminal incidence is not at the path end")
    if word == "?" and terminal is not None and len(graph) == 3 and not adjacency:
        return sp.Integer(1), None, "unique-terminal"
    if "?" in word or adjacency:
        raise UnsupportedReduction("unhandled partial word or condition")
    name = ("S_" if terminal is not None else "C_") + word
    kind = "suffix" if terminal is not None else "motif"
    return sp.Symbol(name), Pattern(name, identity, kind, {"word": word}), kind


def reduce_input(rule) -> Image:
    projected = projected_input(rule)
    if projected is None:
        return Image(sp.Integer(0), {}, {"zero_reason": "contradictory-condition"})
    graph, adjacency, no_edges = projected
    coupled = graph.copy()
    coupled.add_edges_from(no_edges)  # conditions forbid unjustified disconnection
    components, patterns, categories = [], {}, []
    expression = sp.Integer(1)
    pools = Counter()
    for ids in nx.connected_components(coupled):
        ids = set(ids)
        component = graph.subgraph(ids).copy()
        expr, pattern, category = recognize_component(
            component, [a for a in adjacency if a[0] in ids],
            [pair for pair in no_edges if pair[0] in ids and pair[1] in ids])
        if category == "pool":
            pools[pattern.name] += 1
        else:
            expression *= expr
        if pattern:
            patterns[pattern.name] = pattern
        categories.append(category)
        components.append({"vertices": sorted(ids), "factor": str(expr), "justification": category})
    # Disconnected match counts are NOT products in arbitrary graph languages.
    # Free tokens cannot meet incident path positions, tags follow their owners,
    # and gate/feed sorts are disjoint. Refuse cases needing overlap corrections.
    if sum(c in ("motif", "suffix", "unique-terminal") for c in categories) > 1:
        raise UnsupportedReduction("several path components require injective-overlap corrections")
    if "total" in categories and len(components) > 1:
        raise UnsupportedReduction("resource products require factorial/overlap normalization")
    # Distinct free-token components require distinct copies in the host.
    for name, count in pools.items():
        expression *= sp.prod(sp.Symbol(name) - i for i in range(count))
    return Image(sp.expand(expression), patterns, {"component_factorization": components})


class PathCounter:
    def __init__(self, *, include_empty=False, record_compositions=False, state_space=None):
        self.omega = state_space or state_family()
        self.omega.require(state_family(self.omega.definition["alphabet"], self.omega.definition["gated"]))
        self.include_empty = include_empty
        self.record_compositions = record_compositions

    def activity(self, source: SourceRule) -> Image:
        reduced = reduce_input(source.rule)
        reduced.evidence["source_gml"] = source.rule.getGMLString()
        return reduced

    def image(self, pattern: Pattern, source: SourceRule) -> Image:
        calculated = self.compose(pattern.rule, source)
        image = reduce_signed_rules(calculated.composite_coefficients,
                                    calculated.composite_representatives, reduce_input)
        image.evidence.update(composition_orders=calculated.composition_orders,
                              include_empty=self.include_empty)
        if self.record_compositions:
            image.evidence.update(raw_composites=calculated.raw_composites,
                conditional_classes=[{"coefficient": calculated.composite_coefficients.get(key, 0),
                                      "gml": calculated.composite_representatives[key].getGMLString()}
                                     for key in sorted(calculated.composite_representatives)])
        return image

    reduce = staticmethod(reduce_input)

    def compose(self, observable, source):
        return commutator(observable, source.rule, include_empty=self.include_empty,
            rule_rejection=lambda rule: structural_violation(rule) or structural_violation(rule, "right"),
            record_compositions=self.record_compositions)


def derive(motifs=("ABA",), *, alphabet=("A", "B"), gated=False,
           maximum_patterns=64, maximum_expansions=64, include_empty=False,
           record_compositions=False, state_space=None):
    expected = state_family(alphabet, gated)
    omega = expected if state_space is None else state_space.require(expected)
    alphabet, gated = omega.definition["alphabet"], omega.definition["gated"]
    if not motifs or len(set(motifs)) != len(motifs) or any(set(w) - set(alphabet) for w in motifs):
        raise ValueError("distinct nonempty motifs over the selected alphabet required")
    sources = source_rules(alphabet, gated, state_space=omega)
    queries = [query_pattern(word) for word in motifs]
    backend = PathCounter(include_empty=include_empty,
                              record_compositions=record_compositions, state_space=omega)
    support = discover_support(sources, {q.name: q for q in queries},
        counter=backend,
        maximum_patterns=maximum_patterns, maximum_expansions=maximum_expansions)
    return support, sources
