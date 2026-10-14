"""Growing directed A/B paths: Ω, feed/append rules and selected motif counts."""
from __future__ import annotations
from dataclasses import dataclass
import json
import molfoundry as mf
from mfobs.core.pattern_support import Pattern
from mfobs.api import RatedRule
from mfobs.state_space import StateSpace, alphabet_symbols


def state_family(alphabet=("A", "B"), gated=False):
    """Declare Ω for this model, independently of its selected queries."""
    if type(gated) is not bool:
        raise ValueError("gated must be a boolean")
    return StateSpace("terminal-paths", {
        "alphabet": alphabet_symbols(alphabet), "gated": gated,
        "path": "one finite rooted directed path, including the empty path",
        "pools": "arbitrary finite disjoint free tokens, one feed marker per symbol",
        "end": "one terminal link and one End vertex",
        "incidence": {"next": "Next", "terminal": "Terminal", "source": "src", "target": "tgt"},
        "gate": "one GateSource/GateTarget pair, with or without its bound edge" if gated else None,
    })


def graph_rule(left_nodes, left_edges, right_nodes=None, right_edges=None,
               *, name="pattern", adjacency=(), no_edges=()):
    right_nodes = left_nodes if right_nodes is None else right_nodes
    right_edges = left_edges if right_edges is None else right_edges
    common = {i: label for i, label in left_nodes.items() if right_nodes.get(i) == label}
    edges = lambda ee: {(min(i, j), max(i, j), label) for i, j, label in ee}
    le, re = edges(left_edges), edges(right_edges)
    node_text = lambda nn: "\n".join(
        f"node [ id {i} label {json.dumps(label)} ]" for i, label in sorted(nn.items()))
    edge_text = lambda ee: "\n".join(
        f"edge [ source {i} target {j} label {json.dumps(label)} ]" for i, j, label in sorted(ee))
    text = f"rule [ ruleID {json.dumps(name)} context [ {node_text(common)} {edge_text(le & re)} ]"
    text += f" left [ {node_text({i: x for i, x in left_nodes.items() if i not in common})} {edge_text(le-re)} ]"
    text += f" right [ {node_text({i: x for i, x in right_nodes.items() if i not in common})} {edge_text(re-le)} ]"
    for i, j in sorted(set(tuple(sorted(pair)) for pair in no_edges)):
        text += f" constrainNoEdge [ source {i} target {j} ]"
    for i, op, value in sorted(set(adjacency)):
        text += (f' constrainAdj [ id {i} op "{op}" count {value}'
                 ' nodeLabels [ label "Next" label "Terminal" ]'
                 ' edgeLabels [ label "src" label "tgt" ] ]')
    return mf.Rule.fromGMLString(text + " ]", add=False)


@dataclass(frozen=True)
class SourceRule(RatedRule):
    kind: str = ""
    symbol: str | None = None


def source_rules(alphabet=("A", "B"), gated=False, *, state_space=None) -> tuple[SourceRule, ...]:
    """The depth-zero/unrefined typed projection; independent of all queries.

    The schema matches TerminalWordGraphModel: Position/kind/Symbol tags,
    Next and Terminal edge-nodes with src/tgt ports, and the NoEdge gate.
    Atom-to-projection correspondence remains the existing model obligation.
    """
    expected = state_family(alphabet, gated)
    omega = expected if state_space is None else state_space.require(expected)
    alphabet, gated = omega.definition["alphabet"], omega.definition["gated"]
    result = []
    for a in alphabet:
        feed = graph_rule({0: "Feed:" + a}, [],
                          {0: "Feed:" + a, 1: "Position", 2: "Symbol:" + a},
                          [(1, 2, "kind")], name="feed_" + a)
        result.append(SourceRule(rule=feed, rate="beta_" + a, identifier="feed_" + a, kind="feed", symbol=a))
    for a in alphabet:
        nodes = {0: "Position", 1: "Terminal", 2: "End", 3: "Position", 4: "Symbol:" + a}
        left = [(0, 1, "src"), (1, 2, "tgt"), (3, 4, "kind")]
        right_nodes = {i: label for i, label in nodes.items() if i != 1}
        right_nodes.update({5: "Next", 6: "Terminal"})
        right = [(3, 4, "kind"), (0, 5, "src"), (5, 3, "tgt"), (3, 6, "src"), (6, 2, "tgt")]
        if gated:
            nodes.update({7: "GateSource", 8: "GateTarget"})
            right_nodes.update({7: "GateSource", 8: "GateTarget"})
        append = graph_rule(nodes, left, right_nodes, right, name="insert_" + a,
                            adjacency=[(3, "=", 0)], no_edges=[(7, 8)] if gated else [])
        result.append(SourceRule(rule=append, rate="k_" + a, identifier="insert_" + a, kind="append", symbol=a))
    if gated:
        nodes = {0: "GateSource", 1: "GateTarget"}
        edge = [(0, 1, "bound")]
        result.extend((
            SourceRule(rule=graph_rule(nodes, [], nodes, edge,
                       name="donor_bind", no_edges=[(0, 1)]), rate="a", identifier="donor_bind", kind="bind"),
            SourceRule(rule=graph_rule(nodes, edge, nodes, [],
                       name="donor_unbind"), rate="b", identifier="donor_unbind", kind="unbind"),
        ))
    return tuple(result)


def query_pattern(word: str) -> Pattern:
    """Only selected queries use this constructor; support uses composed graphs."""
    if not word:
        raise ValueError("empty motif")
    nodes, edges, last, next_id = {}, [], None, 0
    for a in word:
        position, tag = next_id, next_id + 1
        next_id += 2
        nodes.update({position: "Position", tag: "Symbol:" + a})
        edges.append((position, tag, "kind"))
        if last is not None:
            nodes[next_id] = "Next"
            edges.extend(((last, next_id, "src"), (next_id, position, "tgt")))
            next_id += 1
        last = position
    name = ("N_" if len(word) == 1 else "C_") + word
    return Pattern(name, graph_rule(nodes, edges, name=name,
                   adjacency=[(0, ">", 0)] if len(word) == 1 else []),
                   "motif", {"word": word})
