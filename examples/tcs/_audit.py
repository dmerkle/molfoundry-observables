"""Detailed checks behind the short TCS composition walkthrough.

Requires the supplied MolFoundry Python build and NetworkX. Public materialized composition is checked against the public iterator with
overlap witnesses; this is an API-consistency check, not an independent
implementation of composition.
The host-count witness is evaluated separately with NetworkX monomorphisms.
The script uses the public composition predicate to read the two published
overlap witnesses. It does not implement the global Chem completion calculus
or reproduce the complete table of restricted chemical compositions.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import platform
from typing import Any

import networkx as nx
from networkx.algorithms.isomorphism import (
    GraphMatcher, categorical_edge_match, categorical_node_match,
)

BEFORE_TO_AFTER = {0: 6, 1: 7, 2: 8, 3: 9, 4: 4, 5: 5,
                   6: 0, 7: 1, 8: 2, 9: 3}


def graph_gml(graph: nx.Graph) -> str:
    """Serialize the independent host witness as an ordinary graph, not a rule."""
    lines = ["graph ["]
    lines.extend(
        f'  node [ id {node} label {json.dumps(attrs["label"])} ]'
        for node, attrs in sorted(graph.nodes(data=True))
    )
    lines.extend(
        f'  edge [ source {source} target {target} '
        f'label {json.dumps(attrs["label"])} ]'
        for source, target, attrs in sorted(graph.edges(data=True))
    )
    return "\n".join([*lines, "]", ""])


def rule_record(rule) -> dict[str, Any]:
    """The span and conditions this audit reads, from MolFoundry's public API."""
    from mfobs.core.rules import added_edges, removed_edges

    def side(nodes, edges):
        return {"node": [{"id": v, "label": label} for v, label in nodes],
                "edge": [{"source": u, "target": v, "label": label} for u, v, label in edges]}

    return {
        "left": side([(v.id, v.leftLabel) for v in rule.vertices
                      if v.leftLabel is not None and v.membership.value != 'context'],
                     sorted(removed_edges(rule))),
        "context": side([(v.id, v.leftLabel) for v in rule.vertices if v.membership.value == 'context'],
                        [(e.source.id, e.target.id, e.leftLabel) for e in rule.edges if e.isContext]),
        "right": side([(v.id, v.rightLabel) for v in rule.vertices
                       if v.rightLabel is not None and v.membership.value != 'context'],
                      sorted(added_edges(rule))),
        "matchConstraintsNoEdge": [{"source": c.source, "target": c.target}
                                   for c in rule.matchConstraintsNoEdge],
        "matchConstraintsAdjacency": list(rule.matchConstraintsAdjacency),
        "matchConstraintsLabelAny": list(rule.matchConstraintsLabelAny),
        "matchConstraintsShortestPath": list(rule.matchConstraintsShortestPath),
        "matchCondition": rule.matchCondition,
    }


def plain_span(data: dict[str, Any], rename: dict[int, int] | None = None):
    def rid(value):
        return rename[int(value)] if rename is not None else int(value)
    return {
        side: {
            "nodes": sorted((rid(v["id"]), v["label"]) for v in data[side]["node"]),
            "edges": sorted((*sorted((rid(e["source"]), rid(e["target"]))), e["label"])
                            for e in data[side]["edge"]),
        }
        for side in ("left", "context", "right")
    }


def input_graph(data):
    graph = nx.Graph()
    for side in ("left", "context"):
        for node in data[side]["node"]:
            graph.add_node(int(node["id"]), label=node["label"])
        for edge in data[side]["edge"]:
            graph.add_edge(int(edge["source"]), int(edge["target"]), label=edge["label"])
    return graph


def count_conditioned_input(data, host):
    if data.get("matchCondition") is not None or any(data.get(name) for name in (
        "matchConstraintsAdjacency", "matchConstraintsLabelAny", "matchConstraintsShortestPath"
    )):
        raise ValueError("independent witness checker supports direct NoEdge conditions only")
    matcher = GraphMatcher(
        host, input_graph(data), node_match=categorical_node_match("label", None),
        edge_match=categorical_edge_match("label", None),
    )
    count = 0
    for host_to_pattern in matcher.subgraph_monomorphisms_iter():
        image = {pattern: host_vertex for host_vertex, pattern in host_to_pattern.items()}
        if all(not host.has_edge(image[c["source"]], image[c["target"]])
               for c in data["matchConstraintsNoEdge"]):
            count += 1
    return count


def witness_host(before):
    host = input_graph(before)
    host.add_edge(0, 5, label="-")
    valence, order = {"C": 4, "O": 2, "H": 1}, {"-": 1, "=": 2}
    for node, attrs in list(host.nodes(data=True)):
        deficit = valence[attrs["label"]] - sum(order[e["label"]] for *_, e in host.edges(node, data=True))
        if deficit < 0:
            raise ValueError("overvalent witness atom")
        for _ in range(deficit):
            hydrogen = max(host.nodes) + 1
            host.add_node(hydrogen, label="H")
            host.add_edge(node, hydrogen, label="-")
    assert all(sum(order[e["label"]] for *_, e in host.edges(v, data=True)) == valence[a["label"]]
               for v, a in host.nodes(data=True))
    return host


def verify(input_texts, reaction, queries, compositions) -> dict[str, Any]:
    """Check raw overlaps, multiplicity counters, and one separating host."""
    import molfoundry as mf
    from mfobs.core.composition import witnessed_compositions

    variants = {}
    selected_data = {}
    for name, query in queries.items():
        guarded = name == "conditioned"
        counts, selected, diagnostics = {}, {}, {}
        for side, first, second, overlap in (
            ("before", query, reaction, [[4, 4], [5, 5]]),
            ("after", reaction, query, [[1, 5], [2, 4]]),
        ):
            witnessed = list(witnessed_compositions(first, second, include_empty=False))
            keys = [rule.canonKey for rule, _ in witnessed]
            chosen = [rule for rule, actual in witnessed
                      if sorted(actual) == sorted(tuple(pair) for pair in overlap)]
            if len(chosen) != 1:
                raise AssertionError('selected published overlap is not a single public product')
            raw_result = {'keys': keys, 'candidate_count': None,
                          'selected': {'key': chosen[0].canonKey, 'gml': chosen[0].getGMLString()}}
            public_products = compositions[name][side]
            public_counts = Counter(rule.canonKey for rule in public_products)
            raw_counts = Counter(raw_result["keys"])
            if public_counts != raw_counts:
                raise AssertionError(f"{name}/{side}: public/raw multiplicity counters differ")
            value = raw_result["selected"]
            composite = mf.Rule.fromGMLString(value["gml"], add=False)
            if composite.canonKey != value["key"]:
                raise AssertionError("selected composite GML does not round-trip")
            value["data"] = rule_record(composite)
            if public_counts[value["key"]] != 1:
                raise AssertionError("selected overlap is not a single retained contribution")
            selected[side] = {**value, "overlap_first_to_second": overlap}
            counts[side] = public_counts
            diagnostics[side] = {
                "candidate_overlaps": raw_result["candidate_count"],
                "accepted_raw": sum(raw_counts.values()),
                "canonical_classes": len(raw_counts),
                "public_raw_counter_equal": True,
                "canonical_multiplicity_counter_sha256": sha256(
                    json.dumps(sorted(raw_counts.items()), separators=(",", ":")).encode()
                ).hexdigest(),
            }
        same_positive = plain_span(selected["before"]["data"], BEFORE_TO_AFTER) == plain_span(selected["after"]["data"])
        if not same_positive:
            raise AssertionError("selected positive spans do not match the published correspondence")
        coefficients = {
            side: counts["after"][value["key"]] - counts["before"][value["key"]]
            for side, value in selected.items()
        }
        expected = {"before": -1, "after": 1} if guarded else {"before": 0, "after": 0}
        assert coefficients == expected
        same_conditional = selected["before"]["key"] == selected["after"]["key"]
        assert same_conditional is (not guarded)
        selected_data[name] = selected
        variants[name] = {
            "positive_spans_equal_under_map": same_positive,
            "before_to_after_vertex_map": BEFORE_TO_AFTER,
            "conditional_keys_equal": same_conditional,
            "selected_net_commutator_coefficients": coefficients,
            "enumeration": diagnostics,
            "selected_composites": {
                side: {"gml": value["gml"], "overlap_first_to_second": value["overlap_first_to_second"],
                       "key_sha256": sha256(value["key"].encode()).hexdigest(),
                       "noedge": value["data"]["matchConstraintsNoEdge"]}
                for side, value in selected.items()
            },
        }
    host = witness_host(selected_data["conditioned"]["before"]["data"])
    host_counts = {
        name: {side: count_conditioned_input(value["data"], host)
               for side, value in selected.items()}
        for name, selected in selected_data.items()
    }
    assert host_counts == {"unconditioned": {"before": 1, "after": 1},
                           "conditioned": {"before": 0, "after": 1}}
    engine_js = Path(mf.__file__).parent / "bin" / "molfoundry.js"
    return {
        "schema_version": 2, "status": "PASS",
        "scope": "Local r180 positive-span and transported-condition comparison, with an independent input-match witness. Full global Chem constraints and the full Table F1 computation are not imposed.",
        "article": {"doi": "10.1016/j.tcs.2021.07.026", "long_version": "arXiv:2106.02573",
                    "locations": ["Section 8.3.1", "Section 8.3.2, Figure 6 and r180 discussion", "Section 8.3.3"]},
        "runtime": {"python": platform.python_version(), "platform": platform.platform(),
                    "molfoundry": importlib.metadata.version("molfoundry"),
                    "networkx": nx.__version__, "engine_js_sha256": sha256(engine_js.read_bytes()).hexdigest()},
        "method": "Public RCEvaluator onlyUnique=False counters versus raw overlap counters from the same engine; separately, NetworkX injective non-induced host matching with direct NoEdge evaluation.",
        "composition_convention": "Python first argument executes first; signed commutator is after minus before.",
        "inputs": input_texts,
        "variants": variants,
        "host_witness": {
            "representation": "One valence-complete host graph for the independent match-count witness; not a rewriting rule.",
            "nodes": [{"id": v, "label": a["label"]} for v, a in sorted(host.nodes(data=True))],
            "edges": [{"source": a, "target": b, "label": data["label"]}
                      for a, b, data in sorted(host.edges(data=True))],
            "gml": graph_gml(host),
            "valence_complete": True, "added_distinguishing_edge": [0, 5],
            "input_match_counts": host_counts,
        },
    }
