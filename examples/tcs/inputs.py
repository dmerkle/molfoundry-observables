"""Local covalent restrictions for the aldol rule-composition example."""
from typing import Mapping

from mfobs.core.graph import edge_rows, node_labels
from mfobs.core.gml import rule_side
from mfobs.state_space import covalent_compositions as state_family

def _side(rule, side: str, bond_orders: Mapping[str, int]):
    graph = rule_side(rule, side)
    nodes = node_labels(graph)
    incident = {node: [] for node in nodes}
    valence = {node: 0 for node in nodes}
    for source, target, label in edge_rows(graph):
        order = int(bond_orders[label])
        incident[source].append((label, nodes[target]))
        incident[target].append((label, nodes[source]))
        valence[source] += order
        valence[target] += order
    return nodes, incident, valence


def chemical_rule(rule, chemistry: Mapping) -> bool:
    """The generic check used by the article's ``checkRule`` function."""
    from mfobs.state_space import StateSpace
    if isinstance(chemistry, StateSpace):
        constraints = chemistry.chemical_constraints()
        chemistry = {"bond_orders": constraints["bondOrders"],
                     "atom_labels": constraints["atomLabels"],
                     "max_valence": constraints["maxValence"]}
    bond_orders = chemistry["bond_orders"]
    left = _side(rule, "left", bond_orders)
    right = _side(rule, "right", bond_orders)
    left_nodes, _, left_valence = left
    right_nodes, _, right_valence = right
    if left_nodes != right_nodes or left_valence != right_valence:
        return False
    labels = set(chemistry["atom_labels"])
    maxima = chemistry["max_valence"]
    return all(
        label in labels and left_valence[node] <= int(maxima[label])
        for node, label in left_nodes.items()
    )
