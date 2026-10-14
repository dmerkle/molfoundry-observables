"""Construction and GML serialization of native MolFoundry graphs and conditions."""
from collections import Counter
import json
import networkx as nx
import molfoundry as mf
from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows
from mfobs.core.rules import make_rule, rule_graph
def shift(g,offset):
    return make_graph({u+offset:l for u,l in node_rows(g)},[(u+offset,v+offset,l) for u,v,l in edge_rows(g)])

def shift_rule(rule, offset):
    """Shift every reference in the NoEdge/local-adjacency condition fragment."""
    from mfobs.core.conditions import engine_conditions, engine_no_edges
    left, right = rule_graph(rule, 'left'), rule_graph(rule, 'right')
    mapping = {node: node + offset for node in node_labels(left)}
    return make_rule(rule.name, shift(left, offset), shift(right, offset),
                     tuple((u + offset, v + offset) for u, v in engine_no_edges(rule)),
                     tuple(c.reindex(mapping) for c in engine_conditions(rule)))


def add_rule_context(rule, context):
    """Add a disjoint preserved graph, retaining the supported conditions."""
    from mfobs.core.conditions import engine_conditions, engine_no_edges
    offset = max(node_labels(context), default=-1) + 1 - min(node_labels(rule_graph(rule)), default=0)
    shifted = shift_rule(rule, offset)
    return make_rule(rule.name, union(context, rule_graph(shifted, 'left')),
                     union(context, rule_graph(shifted, 'right')),
                     engine_no_edges(shifted), engine_conditions(shifted))


def union(*gs):
    nodes={};edges=[]
    for g in gs:
        if set(nodes)&set(node_labels(g)): raise ValueError('overlapping IDs')
        nodes.update(node_labels(g));edges.extend(edge_rows(g))
    return make_graph(nodes,edges)

def disjoint_union(gs):
    out=make_graph({}); off=0
    for g in gs:
        # preserve relative input IDs but reserve a disjoint range
        h=shift(g,off-min(node_labels(g),default=0));out=union(out,h);off=max(node_labels(out),default=-1)+1
    return out

def rule_gml(rule: mf.Rule) -> str:
    return rule.getGMLString()


def bond_order_sums(graph: mf.Graph) -> Counter:
    sums = Counter()
    for u,v,label in edge_rows(graph):
        order={'-':1, '=':2}[label]; sums[u]+=order; sums[v]+=order
    return sums

def graph_gml(graph: mf.Graph) -> str:
    """GML retaining external IDs (native getGMLString uses internal IDs)."""
    lines = [f'node [ id {u} label {json.dumps(label)} ]' for u, label in node_rows(graph)]
    lines += [f'edge [ source {u} target {v} label {json.dumps(label)} ]'
              for u, v, label in edge_rows(graph)]
    return 'graph [\n' + '\n'.join(lines) + '\n]\n'


def component_graphs(graph: mf.Graph) -> list[mf.Graph]:
    h=nx.Graph();h.add_nodes_from(node_labels(graph));h.add_edges_from((u,v) for u,v,_ in edge_rows(graph))
    return [make_graph({u:l for u,l in node_rows(graph) if u in ids},
                       [(u,v,l) for u,v,l in edge_rows(graph) if u in ids]) for ids in nx.connected_components(h)]

def nx_graph(graph: mf.Graph) -> nx.Graph:
    """``graph`` as a networkx graph; node and edge labels are the ``label`` attributes."""
    h=nx.Graph()
    h.add_nodes_from((u,{'label':l}) for u,l in node_rows(graph))
    h.add_edges_from((u,v,{'label':l}) for u,v,l in edge_rows(graph))
    return h

def isomorphism_witness(template: mf.Graph, actual: mf.Graph) -> dict[str,int]:
    """A label-preserving isomorphism ``template -> actual``, keyed by ``str`` of the
    template IDs. MolFoundry's ``Graph.isomorphism`` only says whether one exists."""
    matcher=nx.algorithms.isomorphism.GraphMatcher(nx_graph(template),nx_graph(actual),
        node_match=nx.algorithms.isomorphism.categorical_node_match('label',None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match('label',None))
    if not matcher.is_isomorphic():
        raise AssertionError('the graphs are not isomorphic')
    return {str(u):v for u,v in matcher.mapping.items()}

def rule_side(rule, side):
    return rule_graph(rule, side)


def _label_list(key, labels):
    return f' {key} [ '+' '.join('label '+json.dumps(l) for l in labels)+' ]' if labels else ''

def constraint_gml(c):
    """GML of one ``molfoundry.constraints`` match constraint; vertices are named by id."""
    kind=type(c).__name__
    if kind=='NoEdgeConstraint':
        return f'constrainNoEdge [ source {c.source} target {c.target} ]'
    if kind=='AdjacencyConstraint':
        return (f'constrainAdj [ id {c.id} op {json.dumps(c.op.value)} count {c.count}'
                f'{_label_list("nodeLabels",c.nodeLabels)}{_label_list("edgeLabels",c.edgeLabels)} ]')
    if kind=='ShortestPathConstraint':
        return (f'constrainShortestPath [ source {c.source} target {c.target}'
                f' op {json.dumps(c.op.value)} length {c.length} ]')
    if kind=='LabelAnyConstraint':
        return f'constrainLabelAny [ label {json.dumps(c.label)}{_label_list("labels",c.labels) or " labels [ ]"} ]'
    raise TypeError(f'not a MolFoundry match constraint: {kind}')

def condition_gml(c):
    """GML of a ``molfoundry.constraints`` application condition; vertices are named by id."""
    kind=type(c).__name__
    if kind=='ConditionTrue':
        return 'and [ ]'
    if kind in ('ConditionAnd','ConditionOr'):
        return ('and' if kind=='ConditionAnd' else 'or')+' [ '+' '.join(map(condition_gml,c.operands))+' ]'
    if kind=='ConditionNot':
        return f'not [ {condition_gml(c.condition)} ]'
    if kind=='ConditionExists':
        label=lambda x:'' if x.label is None else ' label '+json.dumps(x.label)
        parts=[f'node [ id {v.id}{label(v)} ]' for v in c.vertices]
        parts+=[f'edge [ source {e.source} target {e.target}{label(e)} ]' for e in c.edges]
        if type(c.condition).__name__!='ConditionTrue':
            parts.append(condition_gml(c.condition))
        return 'exists [ '+' '.join(parts)+' ]'
    return constraint_gml(c)

def apply_molfoundry_rule(components, rule, *, graph_database=None):
    dg=mf.DG(graphDatabase=list(graph_database or components))
    with dg.build() as builder:
        return builder.apply(list(components),rule,mode=mf.ApplyMode.UseAny,verbosity=0)


def native_graph_data(graph):
    """Native vertices/edges serialized with their internal vertex IDs."""
    return {"nodes": sorted([[v.id, v.stringLabel] for v in graph.vertices]),
            "edges": sorted([[min(e.source.id, e.target.id),
                              max(e.source.id, e.target.id), e.stringLabel]
                             for e in graph.edges])}


def read_rule_record(gml):
    """Read a native MolFoundry rule, retaining all application conditions."""
    return mf.Rule.fromGMLString(gml, add=False)
