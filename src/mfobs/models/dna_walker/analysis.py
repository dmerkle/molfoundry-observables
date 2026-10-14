"""Linear observable closure of the supplied DNA walker rules and query."""
from mfobs.core.graph import node_rows, edge_rows, node_labels
from mfobs.core.rules import make_rule, rule_graph
from collections import Counter
from fractions import Fraction
from dataclasses import asdict
from mfobs.core.gml import rule_gml
from mfobs.core.engine import commutator
from mfobs.core.linear_closure import close_linear_span, UnsupportedImage
from mfobs.models.dna_walker.inputs import OMEGA, convert_rule, count_rule, identity
from mfobs.models.dna_walker.counting import input_graph, reduce_count


def analyze(rules, query, *, query_name='E', query_scale=1, maximum_dimension=24, state_space=None):
    """Compose literal rules with ONE selected query; discover residuals afresh.

    Parameters are original directed-multigraph JSON definitions, not a saved
    result or a candidate basis. Rule names merely name independent rate channels.
    """
    if state_space is not None:
        state_space.require(OMEGA)
    for definition in rules:
        rule = convert_rule(definition)
        # A local rewrite must preserve the whole scaffold and its motor.
        fixed = lambda g: ({n:l for n,l in node_rows(g) if l != 'edge_leg'},
            {(u,v,b) for u,v,b in edge_rows(g) if node_labels(g)[u] != 'edge_leg' and node_labels(g)[v] != 'edge_leg'})
        if fixed(rule_graph(rule)) != fixed(rule_graph(rule, "right")):
            raise UnsupportedImage('Rules must preserve the scaffold and motor')
        for side in (rule_graph(rule), rule_graph(rule, "right")):
            reduced, _ = reduce_count(side)
            if reduced is None or sum(l == 'chain' for l in node_labels(side).values()) != 2:
                raise UnsupportedImage('Rules must connect valid configurations on one scaffold edge')
    reactions = {r['name']: convert_rule(r) for r in rules}
    if len(reactions) != len(rules):
        raise ValueError('Rule names must be unique')
    representatives, graphs = {}, {}
    def register(graph):
        obj = count_rule(graph)
        key = obj.canonKey
        representatives[key], graphs[key] = obj, graph
        return key
    graph, _ = reduce_count(rule_graph(identity(query, query_name)))
    if graph is None:
        raise ValueError('Selected query is zero on Omega')
    seed = {register(graph): Fraction(query_scale)}
    expressions = {query_name: seed}
    calls = raw = 0
    def image(source, channel):
        nonlocal calls, raw
        vector = Counter()
        trace = []
        for key, scale in expressions[source].items():
            result = commutator(representatives[key], reactions[channel])
            calls += 2
            raw += result.left_raw + result.right_raw
            for term, coefficient in sorted(result.observable_coefficients.items()):
                projected = result.observable_representatives[term]
                full = input_graph(projected)
                reduced, reason = reduce_count(full)
                target = register(reduced) if reduced is not None else None
                if target is not None:
                    vector[target] += scale * coefficient
                trace.append({'coefficient':str(scale * coefficient),
                              'projection':projected.getGMLString(),
                              'input_gml':rule_gml(make_rule('input',full,full)),
                              'reduced_key':target, 'reason':reason})
        return {k:v for k,v in vector.items() if v}, {'reductions':trace}
    closure = close_linear_span(seeds=expressions.copy(), channels=reactions, image=image,
        maximum_dimension=maximum_dimension,
        on_add=lambda name,vector,*_: expressions.__setitem__(name,vector))
    result = asdict(closure)
    result['basis'] = {n:{k:str(c) for k,c in v.items()} for n,v in closure.basis.items()}
    result['patterns'] = {k:rule_gml(make_rule('count',g,g)) for k,g in graphs.items()}
    result['ordered_composition_calls'] = calls
    result['raw_composites'] = raw
    result['synthesis_state_evaluations'] = 0
    return result
