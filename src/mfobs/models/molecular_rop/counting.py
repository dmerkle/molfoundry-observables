"""Exact molecular-pattern reduction, including finite local zero-witness covers."""
from __future__ import annotations
from mfobs.core.graph import edge_rows, has_edge, make_graph, node_labels, node_rows
from mfobs.core.rules import rule_vertex_ids
from collections import Counter
from dataclasses import asdict
from itertools import product
import sympy as sp
from mfobs.core import engine as eng
from mfobs.core.gml import component_graphs, rule_side

from mfobs.core.mixtures import fast_count
from mfobs.core.pattern_support import Image, Pattern, UnsupportedReduction, discover_support, reduce_signed_rules
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.molecular_rop.inputs import _identity, forbidden_specification, query_pattern, source_rules


def _isomorphic(first, second):
    from mfobs.core.molfoundry_execution import canonical_key
    if (first.numVertices, first.numEdges) != (second.numVertices, second.numEdges):
        return False
    return not node_rows(first) or canonical_key(first) == canonical_key(second)


class MolecularReducer:
    def __init__(self, selection, *, maximum_context_depth=8, maximum_candidates=2000000, state_space=None):
        from mfobs.models.molecular_rop.inputs import state_family as molecular_paths
        self.selected = rop.configuration(selection)
        expected = molecular_paths(self.selected.alphabet, self.selected.gated)
        self.omega = expected if state_space is None else state_space.require(expected)
        if not isinstance(maximum_context_depth, int) or isinstance(maximum_context_depth, bool) or maximum_context_depth < 1:
            raise ValueError("local-context depth limit must be a positive integer")
        self.maximum_context_depth = maximum_context_depth
        if not isinstance(maximum_candidates,int) or isinstance(maximum_candidates,bool) or maximum_candidates < 0:
            raise ValueError("candidate overlap limit must be a nonnegative integer")
        self.maximum_candidates = maximum_candidates
        self.specification = forbidden_specification(self.selected, maximum_candidates=maximum_candidates)
        self.component_cache = {}
        self.context_hosts = {}
        self.failure = None

    def _context_hosts(self, bound):
        if bound in self.context_hosts:
            return self.context_hosts[bound]
        # Include all short chains as well as both boundary contexts.
        hosts = [rop.chain(tuple(word)) for length in range(bound + 1)
                 for word in product(self.selected.alphabet, repeat=length)]
        hosts += [rop.monomer(a) for a in self.selected.alphabet]
        hosts += [rop.reservoir(a) for a in self.selected.alphabet]
        if self.selected.gated:
            hosts += [rop.catalyst(False), rop.catalyst(True)]
        for host in hosts:
            self.omega.validate_graph(host)
        self.context_hosts[bound] = tuple(hosts)
        return self.context_hosts[bound]

    def _zero_certificate(self, graph):
        bound = (graph.numVertices - 1) // 6 + 2
        if bound > self.maximum_context_depth:
            raise UnsupportedReduction(f"local-context depth {bound} exceeds limit {self.maximum_context_depth}")
        for host in self._context_hosts(bound):
            if fast_count(graph, host):
                return None
        return {"kind": "no-positive-embedding-in-complete-local-context-cover",
                "vertices": graph.numVertices, "maximum_repeat_span": bound,
                "context_graphs": len(self._context_hosts(bound)),
                "argument": "connected image spans at most floor((vertices-1)/6)+2 repeats; all words up to this length include root/end cases; also test every free-molecule and gate type"}

    def _recognize_positive(self, graph):
        # Recognise a produced graph, rather than seed the discovery queue.
        if graph.numVertices == 1 and node_rows(graph)[0][1] in {"Feed-" + a for a in self.selected.alphabet}:
            return sp.Integer(1), None, "unique-feed-marker"
        oh = make_graph({0: "O", 1: "H"}, [(0, 1, "-")])
        if _isomorphic(graph, oh):
            return sp.Integer(1), None, "unique-terminal-OH"
        for a in self.selected.alphabet:
            if _isomorphic(graph, rop.lactone_core(a)):
                name = "M_" + a
                return sp.Symbol(name), Pattern(name, _identity(name, graph), "pool", {"word": a, "representation": "atom-resolved"}), "free-lactone"
        if self.selected.gated:
            for closed in (False, True):
                if _isomorphic(graph, rop.catalyst(closed)):
                    # Open catalyst positive pattern also matches the closed gate.
                    return sp.Integer(1) if not closed else 1-sp.Symbol("X"), None, "gate-closed" if closed else "gate-positive"
        # Carbonyls determine the number of complete repeats to recognise.
        labels = node_labels(graph)
        n = sum(bond == "=" and {labels[i], labels[j]} == {"C", "O"} for i,j,bond in edge_rows(graph))
        if n and n <= self.maximum_context_depth:
            counts = Counter(labels.values())
            for letters in product(self.selected.alphabet, repeat=n):
                if counts["C"] != sum(1+rop.SPACER_CARBONS[a] for a in letters):
                    continue
                for terminal in (False, True):
                    template = rop.repeat_fragment(letters, terminal=terminal)
                    if template.numVertices != graph.numVertices or not _isomorphic(graph, template):
                        continue
                    word = "".join(letters)
                    name = ("S_" if terminal else "N_" if n == 1 else "C_") + word
                    kind = "suffix" if terminal else "motif"
                    return sp.Symbol(name), Pattern(name, _identity(name, graph), kind, {"word": word, "representation": "atom-resolved"}), "chain-suffix" if terminal else "chain-motif"
        return None

    def _component(self, graph, no_edges):
        identity = _identity("component", graph, no_edges)
        key = identity.canonKey
        if key in self.component_cache:
            return self.component_cache[key]
        if any(has_edge(graph, i, j) for i,j in no_edges):
            result = (sp.Integer(0), {}, {"kind": "contradictory-positive-edge-and-NoEdge"}, "zero")
        else:
            recognized = self._recognize_positive(graph)
            if recognized is None:
                certificate = self._zero_certificate(graph)
                if certificate is None:
                    raise UnsupportedReduction("positive molecular component outside the supported whole-repeat/pool/end grammar")
                result = (sp.Integer(0), {}, certificate, "zero")
            else:
                expression, pattern, category = recognized
                factors = {} if pattern is None else {pattern.name: pattern}
                if category.startswith("gate-"):
                    gate_pairs = [(i,j) for i,j in no_edges if {node_labels(graph)[i],node_labels(graph)[j]} == {"Pd","N"}]
                    if set(no_edges) != set(gate_pairs):
                        raise UnsupportedReduction("unsupported molecular gate absence condition")
                    if category == "gate-positive" and gate_pairs:
                        expression = sp.Symbol("X")
                    if expression.free_symbols:
                        gate = _identity("X", rop.catalyst(False), [(0,4)])
                        factors["X"] = Pattern("X", gate, "gate", {"representation": "atom-resolved"})
                # Recognised organic whole-unit patterns are induced on their
                # images in this constructor language, so other missing edges
                # inside them are necessarily absent. No arbitrary NAC is erased.
                result = (expression, factors, {"kind": category,
                    "NoEdge_reduction": "gate predicate" if category.startswith("gate-") else "induced whole-molecule/unit image in constructor language",
                    "no_edges": list(no_edges), "normalization": "one injective map per oriented occurrence"}, category)
        self.component_cache[key] = result
        return result

    def reduce(self, rule):
        if rule_vertex_ids(rule, "deleted"):
            raise UnsupportedReduction("molecular DPO vertex deletion is not supported; interface cannot be erased")
        if rule.matchCondition is not None or any(getattr(rule, k) for k in ("matchConstraintsAdjacency", "matchConstraintsLabelAny", "matchConstraintsShortestPath")):
            raise UnsupportedReduction("transported molecular condition outside supported NoEdge fragment")
        graph = rule_side(rule, "left")
        no_edges = [tuple(sorted((c.source, c.target))) for c in rule.matchConstraintsNoEdge]
        pieces = component_graphs(graph)
        reductions = []
        # Retain unsupported components until checking whether another factor
        # certifies that the complete pattern is impossible.
        failure = None
        for part in pieces:
            internal = [p for p in no_edges if set(p) <= set(node_labels(part))]
            try:
                value = self._component(part, internal)
            except UnsupportedReduction as exc:
                failure = exc
                continue
            if value[0] == 0:
                return Image(sp.Integer(0), {}, {"reduction": value[2]})
            reductions.append((part, value))
        if failure:
            raise failure
        categories = [v[3] for _,v in reductions]
        if sum(c in ("chain-suffix", "unique-terminal-OH") for c in categories) > 1:
            return Image(sp.Integer(0), {}, {"reduction": {
                "kind": "two-injective-copies-of-the-unique-terminal-OH"}})
        if sum(c.startswith("chain-") or c == "unique-terminal-OH" for c in categories) > 1:
            raise UnsupportedReduction("several chain factors require injective-overlap corrections")
        if sum(c.startswith("gate-") for c in categories) > 1:
            return Image(sp.Integer(0), {}, {"reduction": {
                "kind": "two-injective-copies-of-the-unique-gate"}})
        markers = [node_rows(part)[0][1] for part,value in reductions if value[3] == "unique-feed-marker"]
        if len(set(markers)) != len(markers):
            return Image(sp.Integer(0), {}, {"reduction": {"kind": "two-injective-copies-of-one-unique-marker"}})
        expression, factors, evidence, pools = sp.Integer(1), {}, [], {}
        for part, (value, patterns, argument, category) in reductions:
            if category == "free-lactone":
                name = next(iter(patterns))
                pools[name] = pools.get(name, 0) + 1
            else:
                expression *= value
            factors.update(patterns)
            evidence.append({"vertices": list(node_labels(part)), "factor": str(value), **argument})
        for name, count in pools.items():
            expression *= sp.prod(sp.Symbol(name) - i for i in range(count))
        # Different admitted roles map to different host components (pool,
        # chain, gate, markers); their cross-NoEdge predicates are automatic.
        return Image(sp.expand(expression), factors, {"component_factorization": evidence,
            "factorization_argument": "one chain factor, falling factorials of distinct free-molecule copies, unique disjoint gate/feed components; cross-component NoEdge is automatic"})

    def activity(self, source):
        try:
            reduced = self.reduce(eng.dpo_jump_closure(source.rule))
            reduced.evidence["source_gml"] = source.rule.getGMLString()
            return reduced
        except UnsupportedReduction as exc:
            self.failure = {"stage": "activity", "source_rule": source.identifier,
                            "source_gml": source.rule.getGMLString(), "reason": str(exc)}
            raise

    def image(self, pattern, source):
        calculated = self.compose(pattern.rule, source)
        try:
            image = reduce_signed_rules(calculated.observable_coefficients,
                                        calculated.observable_representatives, self.reduce)
        except UnsupportedReduction as exc:
            term = exc.term
            self.failure = {"stage": "projected-composite-reduction", "pattern": pattern.name,
                "source_rule": source.identifier, "query_gml": pattern.rule.getGMLString(),
                "source_gml": source.rule.getGMLString(), "projected_gml": term["gml"],
                "coefficient": term["coefficient"], "reason": str(exc)}
            raise
        image.evidence.update(
            actual_operands={"observable_gml": pattern.rule.getGMLString(),
                             "reaction_gml": source.rule.getGMLString()},
            representation="atom-resolved molecular DPO", include_empty=True,
            composition_orders={"before": asdict(calculated.left_restricted),
                                "after": asdict(calculated.right_restricted)},
        )
        return image

    def compose(self, observable, source):
        try:
            return eng.commutator(observable, source.rule,
                forbidden_specification=self.specification, include_empty=True)
        except Exception as exc:
            if "COMPOSITION_CANDIDATE_LIMIT" not in str(exc):
                raise
            self.failure = {"stage": "composition", "kind": "candidate-overlap-limit", "pattern": observable.name,
                "source_rule": source.identifier, "query_gml": observable.getGMLString(),
                "source_gml": source.rule.getGMLString(), "maximum_candidates_per_order": self.maximum_candidates,
                "reason": "exact common-overlap composition exceeded its candidate budget; no partial coefficient used"}
            raise UnsupportedReduction(self.failure["reason"]) from exc


def derive(selection, motifs, *, maximum_patterns=64, maximum_expansions=64, maximum_context_depth=8, maximum_candidates=2000000, state_space=None):
    from mfobs.models.molecular_rop.inputs import state_family as molecular_paths
    selected = rop.configuration(selection)
    expected = molecular_paths(selected.alphabet, selected.gated)
    omega = expected if state_space is None else state_space.require(expected)
    if not motifs or len(set(motifs)) != len(motifs) or any(not w or set(w)-set(selected.alphabet) for w in motifs):
        raise ValueError("distinct nonempty query words over the selected alphabet required")
    sources = source_rules(selected, state_space=omega)
    backend = MolecularReducer(selected, maximum_context_depth=maximum_context_depth, maximum_candidates=maximum_candidates, state_space=omega)
    queries = [query_pattern(w) for w in motifs]
    support = discover_support(sources, {p.name:p for p in queries}, counter=backend, maximum_patterns=maximum_patterns, maximum_expansions=maximum_expansions)
    if backend.failure:
        support.trace.append({"minimal_reproducer": backend.failure})
        if backend.failure.get("kind") == "candidate-overlap-limit":
            support.status = "UNKNOWN_COMPOSITION_LIMIT"
    return support, sources
