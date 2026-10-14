"""Reduce transported conditions to exact ring-local Boolean count polynomials."""
from mfobs.core.graph import edge_rows, make_graph, node_labels
from mfobs.core.rules import is_identity, preserves_vertices, rule_graph
from dataclasses import dataclass
from fractions import Fraction
from math import factorial
import sympy as sp
import networkx as nx
from mfobs.core.conditions import OPS, UnsupportedCondition
from mfobs.core.exact_span import ExactIncrementalSpan

from mfobs.core.molfoundry_execution import matches
from mfobs.core.pattern_support import Pattern, Image, UnsupportedReduction, polynomial_record, reduce_signed_rules
from mfobs.core.engine import commutator, dpo_jump_closure
from mfobs.core.polynomial_jump_closure import PolynomialRing
from mfobs.models.chlorination.inputs import query_patterns, _identity


Z = sp.Symbol("Z")


RING = PolynomialRing(tuple(f"h{i+1}" for i in range(4)),
                      boolean_variables=tuple(f"h{i+1}" for i in range(4)))


ZERO, ONE = RING.constant(0), RING.constant(1)


H = tuple(RING.coordinate(name) for name in RING.variables)


@dataclass(frozen=True)
class ConditionalReduction:
    expression: sp.Expr
    evidence: dict


def _symbolic_host(hcl_copies):
    fixed = {0: "C-", **{i: "C" for i in range(1, 6)},
             **{i: "H" for i in (10, 11, 12)}}
    dynamic = {6+i: (i, "H", "Cl") for i in range(4)}
    edges = {(min(u, v), max(u, v)): bond for u, v, bond in
             [(i, (i+1) % 5, ":") for i in range(5)] + [(0, 5, "-")] +
             [(i+1, i+6, "-") for i in range(4)] + [(5, i, "-") for i in (10, 11, 12)]}
    for i in range(hcl_copies):
        h, cl = 20+2*i, 21+2*i
        fixed.update({h: "H", cl: "Cl"})
        edges[h, cl] = "-"
    return fixed, dynamic, edges


def _label_weight(label, vertex, fixed, dynamic):
    if label is None:
        return ONE
    if vertex in fixed:
        return ONE if fixed[vertex] == label else ZERO
    index, positive, negative = dynamic[vertex]
    return H[index] if label == positive else RING.subtract(ONE, H[index]) if label == negative else ZERO


def _input_graph(rule):
    left = rule_graph(rule, "left")
    return node_labels(left), list(edge_rows(left))


def _components(nodes, edges):
    topology = nx.Graph()
    topology.add_nodes_from(nodes)
    topology.add_edges_from((u, v) for u, v, _ in edges)
    return [tuple(ids) for ids in nx.connected_components(topology)]


def _hcl_component_count(nodes, edges):
    return sum(len(ids) == 2 and sorted(nodes[i] for i in ids) == ["Cl", "H"] and
               [bond for u, v, bond in edges if u in ids and v in ids] == ["-"]
               for ids in _components(nodes, edges))


def _weighted_matches(nodes, edges, fixed, dynamic, host_edges, anchors=None):
    """MolFoundry supplies the maps; symbolic labels supply their weights.

    Unique temporary colours fix already-bound vertices in existential
    extensions. All remaining vertices share a colour because their actual
    labels are evaluated in the Boolean polynomial ring.
    """
    anchors = anchors or {}
    colours = {target: f"anchor:{node}" for node, target in anchors.items()}
    host = make_graph({v: colours.get(v, "vertex") for v in fixed.keys() | dynamic.keys()},
                      [(u, v, bond) for (u, v), bond in host_edges.items()])
    pattern = make_graph({v: f"anchor:{v}" if v in anchors else "vertex" for v in nodes},
                         [(u, v, bond) for u, v, bond in edges if bond is not None])
    for mapping in matches(pattern, host):
        if any(host_edges.get(tuple(sorted((mapping[u], mapping[v])))) is None
               for u, v, bond in edges if bond is None):
            continue
        weight = RING.multiply(*(_label_weight(label, mapping[v], fixed, dynamic)
                                 for v, label in nodes.items() if v not in anchors))
        if not weight.is_zero:
            yield mapping, weight


def _condition_value(condition, mapping, current_ids, fixed, dynamic, host_edges):
    if condition is None or type(condition).__name__ == "ConditionTrue":
        return ONE
    kind = type(condition).__name__
    evaluate = lambda child: _condition_value(child, mapping, current_ids, fixed, dynamic, host_edges)
    if kind == "ConditionNot":
        return RING.subtract(ONE, evaluate(condition.condition))
    if kind == "ConditionAnd":
        return RING.multiply(*(evaluate(child) for child in condition.operands))
    if kind == "ConditionOr":
        return RING.subtract(ONE, RING.multiply(*(RING.subtract(ONE, evaluate(child))
                                               for child in condition.operands)))
    if kind != "ConditionExists":
        raise UnsupportedCondition(f"Unsupported generic condition node {kind!r}")
    added = {v.id: v.label for v in condition.vertices}
    if set(added) & current_ids:
        raise UnsupportedCondition("Condition vertex id is already in scope")
    nodes = {**dict.fromkeys(current_ids), **added}
    edges = [(e.source, e.target, e.label) for e in condition.edges]
    if any(u not in nodes or v not in nodes for u, v, _ in edges):
        raise UnsupportedCondition("Transported condition references an unmapped vertex")
    if any(not set(ids) & current_ids for ids in _components(nodes, edges)):
        raise UnsupportedCondition("Unanchored existential is outside the local reducer")
    complement = ONE
    for extension, weight in _weighted_matches(nodes, edges, fixed, dynamic, host_edges, mapping):
        child = _condition_value(condition.condition, extension, set(nodes), fixed, dynamic, host_edges)
        complement = RING.multiply(complement, RING.subtract(ONE, RING.multiply(weight, child)))
    return RING.subtract(ONE, complement)


def _direct_constraints_hold(rule, mapping, fixed, dynamic, edges):
    if any(tuple(sorted((mapping[c.source], mapping[c.target]))) in edges
           for c in rule.matchConstraintsNoEdge):
        return False
    for c in rule.matchConstraintsAdjacency:
        center, count = mapping[c.id], 0
        for (u, v), bond in edges.items():
            if center not in (u, v) or bond not in set(c.edgeLabels or ()):
                continue
            neighbour = v if u == center else u
            if neighbour in dynamic:
                raise UnsupportedCondition("Adjacency test with a symbolic H/Cl neighbour is outside this reducer")
            count += fixed[neighbour] in set(c.nodeLabels or ())
        if c.op.value not in OPS:
            raise UnsupportedCondition(f"Unsupported adjacency comparison {c.op.value!r}")
        if not OPS[c.op.value](count, c.count):
            return False
    if rule.matchConstraintsLabelAny or rule.matchConstraintsShortestPath:
        raise UnsupportedCondition("LabelAny and ShortestPath are outside this reducer")
    return True


def certify_locality(rule):
    """Prove count/condition locality on every admitted tail and pool size.

    Positive maps and every existential extension must stay in the ring/HCl
    core. This is a syntactic sufficient condition, not a sampled-tail test.
    In particular, mere connectedness to the unique C- anchor is insufficient.
    """
    nodes, edges = _input_graph(rule)
    if rule.matchConstraintsLabelAny or rule.matchConstraintsShortestPath:
        raise UnsupportedCondition("LabelAny and ShortestPath are outside the locality certificate")

    def graph_scope(labels, bonds, inherited=()):
        forced = set(inherited) | {u for u, label in labels.items() if label == "C-"}
        for u, v, bond in bonds:
            if bond == ":":
                forced.update(x for x in (u, v) if labels[x] in ("C", "C-"))
            if bond == "-":
                if labels[u] == "C" and labels[v] == "Cl":
                    forced.add(u)
                if labels[v] == "C" and labels[u] == "Cl":
                    forced.add(v)
        loose = [u for u, label in labels.items() if
                 (label == "C" and u not in forced) or
                 (label == "H" and not any(
                     bond == "-" and u in (a, b) and
                     ((b if a == u else a) in forced or labels[b if a == u else a] == "Cl")
                     for a, b, bond in bonds)) or label not in ("C", "C-", "H", "Cl")]
        if loose:
            raise UnsupportedCondition(f"Nonlocal pattern or existential extension at vertices {loose}")
        return forced

    def required_ring(condition, labels):
        """Ring membership implied by truth; OR intersects, AND unions facts.

        Composition can encode a direct aromatic-degree test as nested exists.
        Never infer membership from a negated existential.
        """
        if condition is None:
            return set()
        kind = type(condition).__name__
        if kind in ("ConditionAnd", "ConditionOr"):
            parts = [required_ring(c, labels) for c in condition.operands]
            if not parts:
                return set()
            return set.union(*parts) if kind == "ConditionAnd" else set.intersection(*parts)
        if kind != "ConditionExists":
            return set()
        extended = {**labels, **{v.id: v.label or "" for v in condition.vertices}}
        forced = set()
        for e in condition.edges:
            u, v = e.source, e.target
            if e.label is None:
                continue
            if e.label == ":":
                forced.update(x for x in (u, v) if extended[x] in ("C", "C-"))
            if e.label == "-":
                if extended[u] == "C" and extended[v] == "Cl":
                    forced.add(u)
                if extended[v] == "C" and extended[u] == "Cl":
                    forced.add(v)
        forced |= required_ring(condition.condition, extended)
        return forced & set(labels)

    # Aromatic degree predicates also force a positive input C onto the ring.
    forced = set()
    for c in rule.matchConstraintsAdjacency:
        if c.nodeLabels is None or c.edgeLabels is None:
            raise UnsupportedCondition("Adjacency counting requires explicit node and edge label filters")
        if (nodes[c.id] in ("C", "C-") and set(c.edgeLabels or ()) == {":"} and
                ((c.op.value in ("=", ">=") and c.count >= 1) or
                 (c.op.value == ">" and c.count >= 0))):
            forced.add(c.id)
    forced |= required_ring(rule.matchCondition, nodes)
    forced = graph_scope(nodes, edges, forced)

    def condition_scope(condition, labels, bonds, ring):
        if condition is None:
            return 0
        kind = type(condition).__name__
        if kind == "ConditionTrue":
            return 0
        if kind == "ConditionNot":
            return condition_scope(condition.condition, labels, bonds, ring)
        if kind in ("ConditionAnd", "ConditionOr"):
            return sum(condition_scope(c, labels, bonds, ring) for c in condition.operands)
        if kind != "ConditionExists":
            raise UnsupportedCondition(f"Unsupported generic condition node {kind!r}")
        ids = [v.id for v in condition.vertices]
        added_edges = [(e.source, e.target, e.label) for e in condition.edges]
        reached = set(labels)
        while True:
            previous = set(reached)
            for u, v, _ in added_edges:
                if {u, v} & reached:
                    reached.update((u, v))
            if previous == reached:
                break
        if not set(ids) <= reached:
            raise UnsupportedCondition("Unanchored existential is outside the local reducer")
        if any(v.label is None for v in condition.vertices):
            raise UnsupportedCondition("Wildcard nodes lack an all-tail locality certificate")
        added_edges = [(u, v, "" if bond is None else bond) for u, v, bond in added_edges]
        extended = {**labels, **{v.id: v.label for v in condition.vertices}}
        extended_edges = bonds + added_edges
        extended_ring = graph_scope(extended, extended_edges,
                                    ring | required_ring(condition.condition, extended))
        return 1 + condition_scope(condition.condition, extended, extended_edges, extended_ring)

    extensions = condition_scope(rule.matchCondition, nodes, edges, forced)
    return {"status": "STRUCTURAL_ALL_TAIL_LENGTHS",
            "forced_ring_carbons": sorted(forced), "checked_existentials": extensions,
            "arbitrary_methane_and_hydrogen_pools": True,
            "argument": "Positive maps and all existential witnesses stay in the ring/HCl core; "
                        "edges and tested neighbour profiles there are independent of tail length and CH4/H2 counts."}


def _symbolic_matches(rule, hcl_copies):
    certify_locality(rule)
    nodes, edges = _input_graph(rule)
    fixed, dynamic, host_edges = _symbolic_host(hcl_copies)
    result = []
    for mapping, weight in _weighted_matches(nodes, edges, fixed, dynamic, host_edges):
        if _direct_constraints_hold(rule, mapping, fixed, dynamic, host_edges):
            accepted = RING.multiply(weight, _condition_value(
                rule.matchCondition, mapping, set(nodes), fixed, dynamic, host_edges))
            if not accepted.is_zero:
                result.append((mapping, accepted))
    return result


def _count_on_symbolic_host(rule, hcl_copies):
    return RING.add(*(weight for _, weight in _symbolic_matches(rule, hcl_copies)))


def _mask(powers):
    return sum((bool(power) << i) for i, power in enumerate(powers))


def _boolpoly_record(value):
    return {"terms": [{"monomial": RING.monomial_text(powers), "coefficient": str(c)}
                      for powers, c in sorted(value.terms, key=lambda term: _mask(term[0]))]}


def coordinate_boolpoly(expression, coordinates):
    polynomial = sp.Poly(sp.expand(expression), *(sp.Symbol(n) for n in coordinates), domain=sp.QQ)
    terms = []
    for powers, coefficient in polynomial.terms():
        factors = [coordinates[name] for name, power in zip(coordinates, powers) for _ in range(power)]
        terms.append(RING.scale(RING.multiply(*factors), Fraction(int(coefficient.p), int(coefficient.q))))
    return RING.add(*terms)


def _monomial_pattern(mask, name):
    """Complete the fixed ring scaffold with the monomial's H pendants."""
    nodes = {0: "C-", **{i: "C" for i in range(1, 5)}}
    edges = [(i, (i+1) % 5, ":") for i in range(5)]
    for i in range(4):
        if mask & (1 << i):
            nodes[6+i] = "H"
            edges.append((i+1, 6+i, "-"))
    body = " ".join(f'node [ id {i} label "{label}" ]' for i, label in nodes.items())
    body += " " + " ".join(f'edge [ source {u} target {v} label "{label}" ]' for u, v, label in edges)
    return _identity(name, body)


class LocalConditionalReducer:
    """Reduce transported conditions and generate missing positive counts."""

    def __init__(self, chlorine_total=4, *, state_space=None):
        from mfobs.models.chlorination.inputs import state_family
        expected = state_family(chlorine_total)
        self.omega = expected if state_space is None else state_space.require(expected)
        self.chlorine_total = self.omega.definition["chlorine_total"]
        self.cache = {}
        self.trace = {}
        self.activity_reduction = None
        self.patterns = query_patterns()
        self.coordinates = {name: _count_on_symbolic_host(p.rule, 0) for name, p in self.patterns.items()}
        self.pattern_keys = {p.rule.canonKey for p in self.patterns.values()}
        self.span = ExactIncrementalSpan(self._vector(ONE))
        for count in self.coordinates.values():
            self.span.add(self._vector(count))

    @staticmethod
    def _vector(value):
        coefficients = {_mask(powers): c for powers, c in value.terms}
        return tuple(coefficients.get(i, Fraction()) for i in range(1 << len(H)))

    def express(self, local, *, discover=True):
        def solve():
            coefficients = self.span.solve(self._vector(local))
            return None if coefficients is None else sp.expand(sp.Rational(coefficients[0]) + sum(
                sp.Rational(c)*sp.Symbol(name) for name, c in zip(self.coordinates, coefficients[1:])))

        expression = solve()
        if expression is not None:
            return expression
        if not discover:
            raise UnsupportedCondition("Joint activities require additional positive patterns")
        for powers, _ in sorted(local.terms, key=lambda term: _mask(term[0])):
            mask = _mask(powers)
            if not mask:
                continue
            name = f"P{len(self.patterns)}"
            rule = _monomial_pattern(mask, name)
            key = rule.canonKey
            if key in self.pattern_keys:
                continue
            count = _count_on_symbolic_host(rule, 0)
            vector = self._vector(count)
            if self.span.solve(vector) is not None:
                continue
            self.span.add(vector)
            self.pattern_keys.add(key)
            self.coordinates[name] = count
            self.patterns[name] = Pattern(name, rule, "positive", {
                "construction": "ring scaffold with H pendants from a conditional-count monomial",
                "gml": rule.getGMLString(), "raw_count": _boolpoly_record(count)})
            expression = solve()
            if expression is not None:
                return expression
        raise UnsupportedCondition("Conditional count is outside the generated positive-pattern span")

    def equivalent(self, left, right):
        return coordinate_boolpoly(left-right, self.coordinates).is_zero

    def reduce(self, rule):
        key = rule.canonKey
        if key in self.cache:
            return self.cache[key]
        if not preserves_vertices(rule):
            raise UnsupportedCondition("Vertex-deleting jumps need dangling-aware counting")
        if not is_identity(rule):
            raise UnsupportedCondition("Reduction requires a jump-projected rule")
        nodes, edges = node_labels(rule_graph(rule, "left")), list(edge_rows(rule_graph(rule, "left")))
        for ids in _components(nodes, edges):
            if any(nodes[i] in ("C", "C-") for i in ids):
                continue
            if _hcl_component_count({i: nodes[i] for i in ids}, [e for e in edges if e[0] in ids]) != 1:
                raise UnsupportedCondition("Only ring-anchored patterns and HCl components are supported")
        copies = _hcl_component_count(nodes, edges)
        if copies > 1:
            raise UnsupportedCondition("This prototype supports at most one explicit HCl component")
        local = RING.scale(_count_on_symbolic_host(rule, copies), Fraction(1, factorial(copies)))
        local_expression = self.express(local)
        multiplicity = sp.prod(self.chlorine_total-Z-offset for offset in range(copies))
        expression = sp.expand(multiplicity*local_expression)
        result = ConditionalReduction(expression, {
            "method": "raw-MolFoundry-condition-AST -> symbolic local matches -> finite positive-pattern expansion",
            "state_space_enumeration": False, "hcl_components": copies,
            "locality": certify_locality(rule), "local_boolean_polynomial": _boolpoly_record(local),
            "positive_expansion": str(local_expression), "hcl_falling_factor": str(sp.expand(multiplicity)),
            "coordinate_polynomial": str(expression), "has_generic_matchCondition": rule.matchCondition is not None})
        self.cache[key] = result
        return result

    def reduce_image(self, rule):
        try:
            reduced = self.reduce(rule)
            return Image(reduced.expression, _factors(reduced.expression, self.patterns), reduced.evidence)
        except UnsupportedCondition as exc:
            raise UnsupportedReduction(str(exc)) from exc

    def activity(self, reaction):
        activity_rule = dpo_jump_closure(reaction.rule, name="context_nac_activity")
        self.activity_reduction = self.reduce_image(activity_rule)
        expr = self.activity_reduction.expression
        return Image(expr, _factors(expr, self.patterns), {
            "method": "DPO activity identity -> raw condition AST reduction",
            "state_space_enumeration": False,
            "reduction": self.activity_reduction.evidence,
        })

    def image(self, pattern, reaction):
        result = commutator(pattern.rule, reaction.rule,
                            connected=False, include_empty=True)
        image = reduce_signed_rules(result.observable_coefficients, result.observable_representatives, self.reduce_image)
        expr = image.expression
        terms = [{"canonical_key": t["key"], "coefficient": t["coefficient"], "reduction": t["reduction"]}
                 for t in image.evidence["terms"]]
        self.trace[pattern.name] = {
            "left_raw_compositions": result.left_raw,
            "right_raw_compositions": result.right_raw,
            "projected_conditional_terms": len(terms),
            "terms": terms,
            "result": polynomial_record(expr),
        }
        return Image(expr, _factors(expr, self.patterns), {
            "method": "raw MolFoundry commutator -> jump terms -> AST reducer",
            "state_space_enumeration": False,
            "left_raw_compositions": result.left_raw,
            "right_raw_compositions": result.right_raw,
            "projected_conditional_terms": len(terms),
        })


def _factors(expression, patterns):
    return {
        str(symbol): patterns[str(symbol)]
        for symbol in expression.free_symbols
        if str(symbol) in patterns
    }
