"""DNA walker Ω and conversion of literal directed scaffold inputs to native rules."""

from mfobs.core.graph import make_graph
from mfobs.core.rules import make_rule
from mfobs.core.conditions import identity_rule
from mfobs.core.linear_closure import UnsupportedImage
from mfobs.state_space import StateSpace


def state_family():
    """Declare Ω for this model, independently of its selected queries."""
    return StateSpace("dna-walker", {
        "motors": 1, "legs": 2, "minimum_cycle_length": 5,
        "scaffold": "one finite directed simple cycle",
        "leg_targets": "coincident or adjacent sites; the two legs remain distinct",
        "incidence": {"source": "src", "target": "tgt", "scaffold": "edge_dna", "leg": "edge_leg"},
    })


OMEGA = state_family()


def convert_rule(rule):
    """Encode the exact, unconditioned directed-multigraph interface in JSON."""
    if set(rule) != {'name', 'left', 'right', 'preserved_nodes', 'preserved_edges'}:
        raise UnsupportedImage('Unsupported rule fields or application conditions')
    for graph in (rule['left'], rule['right']):
        if set(graph) != {'nodes', 'edges'}:
            raise UnsupportedImage('Unsupported graph fields')
        if any(set(n) != {'id', 'label'} or n['label'] not in {'bimotor', 'chain'}
               for n in graph['nodes']):
            raise UnsupportedImage('Expected labelled motor/scaffold object vertices')
        if any(set(e) != {'id', 'source', 'target'} for e in graph['edges']):
            raise UnsupportedImage('Expected unlabelled directed multigraph edges')
        for collection in ('nodes', 'edges'):
            ids = [item['id'] for item in graph[collection]]
            if len(set(ids)) != len(ids):
                raise ValueError('Vertex and edge IDs must each be unique')
    node_map = dict(rule['preserved_nodes'])
    edge_map = dict(rule['preserved_edges'])
    reverse_nodes = {v:k for k,v in node_map.items()}
    reverse_edges = {v:k for k,v in edge_map.items()}
    for kind, forward, reverse in [('nodes', node_map, reverse_nodes), ('edges', edge_map, reverse_edges)]:
        if len(forward) != len(reverse) or len(forward) != len(rule['preserved_'+kind]):
            raise ValueError('Noninjective interface')
        if (set(forward) - {n['id'] for n in rule['left'][kind]}
                or set(reverse) - {n['id'] for n in rule['right'][kind]}):
            raise ValueError('Interface references an absent vertex or edge')
    ids = {}
    def allocate(key):
        if key not in ids:
            ids[key] = len(ids)
        return ids[key]
    def encode(graph, side):
        labels = {n['id']:n['label'] for n in graph['nodes']}
        def key(kind, identifier, reverse):
            if side == 'right' and identifier in reverse:
                return kind, 'left', reverse[identifier]
            return kind, side, identifier
        nodes = {allocate(key('node', n, reverse_nodes)):label for n,label in labels.items()}
        edges = []
        for edge in graph['edges']:
            vertex = allocate(key('edge', edge['id'], reverse_edges))
            nodes[vertex] = 'edge_leg' if labels[edge['source']] == 'bimotor' else 'edge_dna'
            source = allocate(key('node', edge['source'], reverse_nodes))
            target = allocate(key('node', edge['target'], reverse_nodes))
            edges.extend([(source, vertex, 'src'), (vertex, target, 'tgt')])
        return make_graph(nodes, edges)
    return make_rule(rule['name'], encode(rule['left'], 'left'), encode(rule['right'], 'right'))


def identity(graph, name):
    return convert_rule({'name':name, 'left':graph, 'right':graph,
        'preserved_nodes':[[n['id'], n['id']] for n in graph['nodes']],
        'preserved_edges':[[e['id'], e['id']] for e in graph['edges']]})


def count_rule(graph, name='count'):
    return identity_rule(graph, name)
