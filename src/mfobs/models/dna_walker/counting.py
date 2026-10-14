"""Scaffold count reduction and DPO input admissibility for DNA walking."""
from mfobs.core.graph import node_rows, make_graph, edge_rows, node_labels
from collections import Counter

from mfobs.core.conditions import engine_has_conditions
from mfobs.core.linear_closure import UnsupportedImage
from mfobs.core.rules import rule_graph, rule_vertex_ids, removed_edges
from mfobs.models.dna_walker.inputs import OMEGA


MINIMUM_CYCLE_LENGTH = OMEGA.definition["minimum_cycle_length"]


def input_graph(projected):
    """Certify DPO dangling before replacing a jump projection by an input count."""
    if engine_has_conditions(projected):
        raise UnsupportedImage('Additional application conditions are unsupported')
    left = rule_graph(projected)
    nodes, edges = node_labels(left), edge_rows(left)
    deleted = set(rule_vertex_ids(projected, "deleted"))
    for node in deleted:
        ports = [(v if u == node else u, label) for u,v,label in edges if node in (u,v)]
        expected = {'src': 'bimotor' if nodes[node] == 'edge_leg' else 'chain', 'tgt': 'chain'}
        if (nodes[node] not in ('edge_leg', 'edge_dna') or len(ports) != 2
                or sorted(label for _,label in ports) != ['src','tgt']
                or any(nodes[other] != expected[label] for other,label in ports)):
            raise UnsupportedImage('DPO dangling requires complete typed edge vertices')
    if any(u not in deleted and v not in deleted for u, v, _ in removed_edges(projected)):
        raise UnsupportedImage('Only complete edge vertices may be deleted')
    return make_graph(nodes, edges)


def reduce_count(graph):
    """Return an equal-count subpattern, or zero, with an explicit reason.

    On Omega, a connected scaffold path with fewer than five sites extends
    uniquely from its leg anchors. Stripping its unanchored ends is therefore
    a bijection of injective matches, including both permutations of parallel
    legs. Longer/disconnected/incomplete contexts are NOT certified by sampling.
    """
    labels = node_labels(graph)
    if set(labels.values()) - {'bimotor','chain','edge_leg','edge_dna'}:
        raise UnsupportedImage('Unknown incidence type')
    motors = [n for n,l in node_rows(graph) if l == 'bimotor']
    legs = [n for n,l in node_rows(graph) if l == 'edge_leg']
    sites = {n for n,l in node_rows(graph) if l == 'chain'}
    if len(motors) > 1 or len(legs) > 2:
        return None, 'zero: exceeds the motor/leg inventory'
    if len(motors) != 1 or len(legs) != 2:
        raise UnsupportedImage('Reduction requires one motor and both legs')
    for u,v,bond in edge_rows(graph):
        if bond not in ('src','tgt') or ((labels[u].startswith('edge_')) == labels[v].startswith('edge_')):
            raise UnsupportedImage('Not a typed incidence graph')
    decoded = {}
    for node,label in node_rows(graph):
        if not label.startswith('edge_'):
            continue
        ports = [(v if u == node else u,bond) for u,v,bond in edge_rows(graph) if node in (u,v)]
        if any(sum(b == port for _,b in ports) > 1 for port in ('src','tgt')):
            return None, 'zero: an edge has multiple sources or targets'
        if len(ports) != 2:
            raise UnsupportedImage('Incomplete edge incidence')
        endpoints = dict((bond,other) for other,bond in ports)
        if (labels[endpoints['src']] != ('bimotor' if label == 'edge_leg' else 'chain')
                or labels[endpoints['tgt']] != 'chain'):
            return None, 'zero: incompatible endpoint types'
        decoded[node] = (endpoints['src'],endpoints['tgt'])
    dna = {n:pair for n,pair in decoded.items() if labels[n] == 'edge_dna'}
    if any(c > 1 for c in Counter(a for a,b in dna.values()).values()) or any(c > 1 for c in Counter(b for a,b in dna.values()).values()):
        return None, 'zero: scaffold branches or parallel scaffold edges'
    successor = dict(dna.values())
    starts = sites - {b for a,b in dna.values()}
    if not starts and len(sites) < MINIMUM_CYCLE_LENGTH:
        return None, 'zero: scaffold cycle shorter than five'
    if len(sites) >= MINIMUM_CYCLE_LENGTH or len(starts) != 1:
        raise UnsupportedImage('Only connected scaffold paths shorter than five are certified')
    path = []
    node = next(iter(starts))
    while node not in path:
        path.append(node)
        if node not in successor:
            break
        node = successor[node]
    if set(path) != sites or path[-1] in successor:
        raise UnsupportedImage('Disconnected or cyclic scaffold context')
    targets = {decoded[n][1] for n in legs}
    positions = sorted(path.index(n) for n in targets)
    if positions[-1] - positions[0] > 1:
        return None, 'zero: legs are neither coincident nor adjacent'
    keep = set(motors + legs) | targets
    keep |= {n for n,(a,b) in dna.items() if a in targets and b in targets}
    reduced = make_graph({n:labels[n] for n in keep},
                         [(u,v,b) for u,v,b in edge_rows(graph) if u in keep and v in keep])
    return reduced, 'unique extension of anchored scaffold path'
