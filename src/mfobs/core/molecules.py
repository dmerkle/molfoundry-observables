"""Shared construction helpers for explicit covalent graphs.

These functions do not match or apply rules. Hydrogen completion is only for
full molecular hosts; rule patterns retain exactly their declared atoms.
"""
from collections import Counter

from mfobs.core.gml import bond_order_sums, component_graphs
import molfoundry as mf
from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows

VALENCE = {'C': 4, 'N': 3, 'S': 2, 'O': 2, 'H': 1}
BOND_ORDERS = {'-': 1, '=': 2}


def complete_with_h(graph: mf.Graph, *, valence_overrides=None) -> mf.Graph:
    """Complete FULL molecular hosts only; never saturate a rule pattern."""
    overrides = valence_overrides or {}
    nodes, edges = node_labels(graph), list(edge_rows(graph))
    sums = bond_order_sums(graph)
    nxt = max(nodes, default=-1) + 1
    for u, label in node_rows(graph):
        deficit = overrides.get(u, VALENCE[label]) - sums[u]
        if deficit < 0:
            raise ValueError(f'overvalent {label} at {u}')
        for _ in range(deficit):
            nodes[nxt] = 'H'
            edges.append((u, nxt, '-'))
            nxt += 1
    return make_graph(nodes, edges)


def formula(graph: mf.Graph) -> dict[str, int]:
    return dict(sorted(Counter(label for _, label in node_rows(graph)).items()))


def automorphism_count(graph: mf.Graph) -> int:
    """Use MolFoundry's public group API, retaining every automorphism."""

    return len(graph.aut().generate())


def smiles(graph: mf.Graph) -> str:
    """Join native component SMILES, retaining every copy in a molecular mixture."""
    return '.'.join(sorted(part.smiles for part in component_graphs(graph)))
