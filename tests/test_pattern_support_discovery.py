"""Executable discovery and independent finite-host semantics regressions."""
from __future__ import annotations

import mfobs.models.aba.events as d_events
import mfobs.models.aba.inputs as d_inputs
from mfobs.checks import terminal_paths as word_check
from dataclasses import replace
from itertools import product
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import networkx as nx
import sympy as sp

from mfobs.core.pattern_support import Image, Pattern, discover_support, UnsupportedReduction
from mfobs.core.rules import rule_graph, rule_vertex_ids, added_edges, removed_edges
from mfobs.core.graph import node_labels
from mfobs.models.aba import counting as d


def host(word, pools, opened):
    graph = nx.Graph()
    def node(label):
        i = len(graph)
        graph.add_node(i, label=label)
        return i
    def tagged(symbol):
        i = node("Position"); j = node(symbol)
        graph.add_edge(i, j, label="kind")
        return i
    last = tagged("Root")
    for letter in word:
        new = tagged("Symbol:" + letter); e = node("Next")
        graph.add_edge(last, e, label="src");graph.add_edge(e, new, label="tgt")
        last = new
    end = node("End");e = node("Terminal")
    graph.add_edge(last, e, label="src");graph.add_edge(e, end, label="tgt")
    for a, number in pools.items():
        node("Feed:" + a)
        for _ in range(number):tagged("Symbol:" + a)
    u, v = node("GateSource"), node("GateTarget")
    if not opened:graph.add_edge(u, v, label="bound")
    return graph


def admissible_maps(rule, graph):
    pattern = d.side_graph(rule)
    match = lambda a, b: a["label"] == b["label"]
    for inverse in nx.algorithms.isomorphism.GraphMatcher(
        graph, pattern, node_match=match, edge_match=match).subgraph_monomorphisms_iter():
        mapping = {j: i for i, j in inverse.items()}
        if any(graph.has_edge(mapping[c.source], mapping[c.target])
               for c in rule.matchConstraintsNoEdge):continue
        okay = True
        for condition in rule.matchConstraintsAdjacency:
            center = mapping[condition.id]
            number = sum(graph.nodes[v]["label"] in condition.nodeLabels
                         and graph[center][v]["label"] in condition.edgeLabels
                         for v in graph[center])
            if not ((number == condition.count) if condition.op.value == "=" else (number > condition.count)):
                okay = False;break
        if okay:yield mapping


def target(rule, graph, mapping):
    result = graph.copy();mapping = dict(mapping)
    right_labels = node_labels(rule_graph(rule, "right"))
    for u, v, _ in removed_edges(rule):
        result.remove_edge(mapping[u], mapping[v])
    for node in rule_vertex_ids(rule, "deleted"):
        assert result.degree(mapping[node]) == 0, "dangling condition"
        result.remove_node(mapping[node])
    for node in rule_vertex_ids(rule, "changed"):
        result.nodes[mapping[node]]["label"] = right_labels[node]
    next_id = max(graph, default=-1) + 1
    for node in rule_vertex_ids(rule, "created"):
        mapping[node] = next_id
        result.add_node(next_id, label=right_labels[node]);next_id += 1
    for u, v, label in added_edges(rule):
        result.add_edge(mapping[u], mapping[v], label=label)
    return result


class DiscoveryTests(unittest.TestCase):
    def test_generic_control_has_no_candidate_basis(self):
        x,y=sp.symbols("x y");source=SimpleNamespace(identifier="r")
        def image(pattern, rule):
            return Image(x*x+y,{"y":Pattern("y",object(),"arbitrary")}) if pattern.name=="x" else Image(sp.Integer(0),{})
        result=discover_support([source],{"x":Pattern("x",object(),"arbitrary")},
             counter=SimpleNamespace(activity=lambda r:Image(sp.Integer(1),{}), image=image))
        self.assertEqual(result.status,"POLYNOMIAL_SUPPORT_CLOSED")
        self.assertEqual(set(result.patterns),{"x","y"})

    def test_single_aba_query(self):
        with patch.object(d,"query_pattern",wraps=d_inputs.query_pattern) as selected:
            result,rules=d.derive(("ABA",))
            self.assertEqual(selected.call_count,1)
        self.assertEqual(result.status,"POLYNOMIAL_SUPPORT_CLOSED")
        self.assertEqual(set(result.patterns),{"C_ABA","M_A","M_B","S_AB","S_A"})
        self.assertEqual(result.origins["S_AB"]["pattern"],"C_ABA")
        self.assertEqual(result.origins["S_A"]["pattern"],"S_AB")
        self.assertEqual(len(d_events.compile_uniform_events(result,rules)["events"]),8)

    def test_depth_zero_projection_uses_existing_source_schema(self):
        # Literal output of the inspected TerminalWordGraphModel construction:
        # root-free depth-zero terminal fragment, free adjacency, optional gate.
        for gated in (False, True):
            gate_nodes = 'node [ id 5 label "GateSource" ] node [ id 6 label "GateTarget" ]' if gated else ''
            start = 7 if gated else 5
            condition = 'constrainNoEdge [ source 5 target 6 ]' if gated else ''
            gml = f'''rule [ context [ node [ id 0 label "Position" ]
                node [ id 1 label "End" ] node [ id 3 label "Position" ]
                node [ id 4 label "Symbol:A" ] edge [ source 3 target 4 label "kind" ] {gate_nodes} ]
                left [ node [ id 2 label "Terminal" ] edge [ source 0 target 2 label "src" ]
                       edge [ source 2 target 1 label "tgt" ] ]
                right [ node [ id {start} label "Next" ] node [ id {start+1} label "Terminal" ]
                        edge [ source 0 target {start} label "src" ] edge [ source {start} target 3 label "tgt" ]
                        edge [ source 3 target {start+1} label "src" ] edge [ source {start+1} target 1 label "tgt" ] ]
                constrainAdj [ id 3 op "=" count 0 nodeLabels [ label "Next" label "Terminal" ]
                    edgeLabels [ label "src" label "tgt" ] ] {condition} ]'''
            reference = d.mf.Rule.fromGMLString(gml, add=False)
            generated = next(r.rule for r in d_inputs.source_rules(gated=gated) if r.identifier == "insert_A")
            self.assertEqual(reference.canonKey, generated.canonKey)

    def test_gated_aba(self):
        result,rules=d.derive(("ABA",),gated=True)
        self.assertEqual(set(result.patterns),{"C_ABA","M_A","M_B","X","S_AB","S_A"})
        self.assertEqual(str(result.images["C_ABA"]["insert_A"]),"M_A*S_AB*X")
        self.assertEqual(len(d_events.compile_uniform_events(result,rules)["events"]),10)
        self.assertEqual(word_check.validate_words(result,rules,maximum_length=5)["status"],"PASS")

    def test_unseen_motifs_and_renaming(self):
        for word,alphabet in [("ABB",("A","B")),("ABBA",("A","B")),("AAAA",("A","B")),("XYX",("X","Y"))]:
            result,rules=d.derive((word,),alphabet=alphabet)
            self.assertEqual(result.status,"POLYNOMIAL_SUPPORT_CLOSED",result.reason)
            self.assertEqual(word_check.validate_words(result,rules,maximum_length=4)["status"],"PASS")
            self.assertTrue(all(p.kind!="suffix" or p.metadata["word"] in word[:-1]
                                for p in result.patterns.values()))

    def test_all_old_preset_queries(self):
        result,rules=d.derive(("A","B","AA","AB","BA","BB","ABA"))
        self.assertEqual(set(n for n,p in result.patterns.items() if p.kind=="suffix"),{"S_A","S_B","S_AB"})
        self.assertEqual(word_check.validate_words(result,rules,maximum_length=4)["status"],"PASS")

    def test_empty_overlap_cancellation(self):
        left,_=d.derive(("AAA",),include_empty=False)
        right,_=d.derive(("AAA",),include_empty=True)
        self.assertEqual(left.images,right.images)
        self.assertEqual(left.activities,right.activities)

    def test_pattern_budget_retains_obligation(self):
        result,_=d.derive(("ABA",),maximum_patterns=3)
        self.assertEqual(result.status,"UNKNOWN_PATTERN_LIMIT")
        self.assertIn("C_ABA",result.pending)

    def test_expansion_budget(self):
        result,rules=d.derive(("ABA",),maximum_expansions=0)
        self.assertEqual(result.status,"UNKNOWN_EXPANSION_LIMIT")
        self.assertIn("C_ABA",result.pending)
        self.assertEqual(d_events.compile_uniform_events(result,rules)["status"],"UNKNOWN_SUPPORT")

    def test_unsaturated_dpo_deletion_is_rejected(self):
        rule=d_inputs.graph_rule({0:"Position"},[],{},[],name="delete-vertex")
        with self.assertRaises(UnsupportedReduction):d.projected_input(rule)

    def test_factorization_does_not_assume_arbitrary_independence(self):
        n={0:"Position",1:"Symbol:A",2:"Position",3:"Symbol:A"}
        rule=d_inputs.graph_rule(n,[(0,1,"kind"),(2,3,"kind")],adjacency=[(0,"=",0),(2,"=",0)])
        self.assertEqual(d.reduce_input(rule).expression, sp.Symbol("M_A")**2-sp.Symbol("M_A"))

    def test_changed_source_needs_new_uniformity_argument(self):
        result,rules=d.derive(("ABA",))
        changed=list(rules)
        changed[0]=replace(changed[0],rule=d_inputs.graph_rule({0:"Feed:A"},[],name="identity-feed"))
        self.assertEqual(d_events.compile_uniform_events(result,changed)["status"],"UNKNOWN_UNIFORM_MATCH_OBLIGATION")

    def test_literal_graph_counts_and_targets(self):
        result,rules=d.derive(("ABA",),gated=True)
        comparison=lambda a,b:nx.is_isomorphic(a,b,
            node_match=lambda x,y:x["label"]==y["label"],edge_match=lambda x,y:x["label"]==y["label"])
        cases=0
        for length in range(5):
            for letters in product("AB",repeat=length):
                word="".join(letters);pools={"A":2,"B":1}
                for opened in (False,True):
                    graph=host(word,pools,opened)
                    values=word_check.coordinate_values(result.patterns,word,pools,opened)
                    for name,pattern in result.patterns.items():
                        self.assertEqual(sum(1 for _ in admissible_maps(pattern.rule,graph)),values[name])
                    for source in rules:
                        expected_pool=dict(pools);expected_word=word;expected_open=opened
                        if source.kind=="feed":
                            expected_pool[source.symbol]+=1;count=1
                        elif source.kind=="append":
                            count=pools[source.symbol]*int(opened)
                            expected_pool[source.symbol]-=1;expected_word+=source.symbol
                        elif source.kind=="bind":count=int(opened);expected_open=False
                        else:count=int(not opened);expected_open=True
                        mappings=list(admissible_maps(source.rule,graph))
                        self.assertEqual(len(mappings),count)
                        for mapping in mappings:
                            self.assertTrue(comparison(target(source.rule,graph,mapping),host(expected_word,expected_pool,expected_open)))
                        cases+=1
        self.assertEqual(cases,372)


if __name__=="__main__":unittest.main()
