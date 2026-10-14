"""State-language checks for directed graphs in an incidence encoding.

An edge-node joins source and target object-nodes through differently labelled
bonds. MolFoundry performs DPO composition and condition transport. Sorts, port
bounds, endpoint and relation typing and the terminal-successor ban are local
to one vertex: they are passed as MolFoundry incidence constraints, which the
engine checks and uses to prune the overlap search. Only parallel directed
edges and directed cycles need a Python predicate.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import partial
from typing import Mapping
from mfobs.core.composition import CompositionDiagnostics, group_compositions, rule_labels
from mfobs.core.graph import node_labels, edge_rows

IncidenceCompositionDiagnostics = CompositionDiagnostics


class TypedIncidencePredicate:
    """The incidence structural theory as a Python check of one graph.

    ``violation`` is the complete check, the reference semantics of
    ``typed_incidence_constraints`` together with ``directed_violation``; the
    latter is the part that is not local to one vertex.
    """

    def __init__(self, specification: Mapping):
        self.spec = specification
        self.objects = set(specification["objectLabels"])
        self.edge_nodes = set(specification["edgeNodeLabels"])
        self.relations = {bond: {tuple(sorted(pair)) for pair in pairs}
                          for bond, pairs in specification.get("relationBonds", {}).items()}

    def _ports(self, graph):
        """Source and target objects of each edge-node, by incidence bond."""
        labels, spec = node_labels(graph), self.spec
        ports = {u: ([], []) for u, label in labels.items() if label in self.edge_nodes}
        for u, v, bond in edge_rows(graph):
            if bond in (spec["sourceBond"], spec["targetBond"]):
                edge, other = (u, v) if u in ports else (v, u)
                if edge in ports and other not in ports:
                    ports[edge][bond == spec["targetBond"]].append(other)
        return ports

    def directed_violation(self, graph):
        spec = self.spec
        arcs = [(sources[0], targets[0]) for sources, targets in self._ports(graph).values()
                if sources and targets]
        if spec["simpleDirected"]:
            seen = set()
            for source, target in arcs:
                key = (min(source, target), max(source, target)) if spec["simpleUnderlying"] else (source, target)
                if key in seen:
                    return "parallel-directed-edge"
                seen.add(key)
        if spec["acyclicDirected"]:
            successors = defaultdict(list)
            for source, target in arcs:
                successors[source].append(target)
            colour = {}
            for start in list(successors):
                if start in colour:
                    continue
                colour[start] = 1
                stack = [(start, iter(successors[start]))]
                while stack:
                    node, pending = stack[-1]
                    following = next(pending, None)
                    if following is None:
                        colour[node] = 2
                        stack.pop()
                    elif colour.get(following) == 1:
                        return "directed-cycle"
                    elif following not in colour:
                        colour[following] = 1
                        stack.append((following, iter(successors[following])))
        return None

    def violation(self, graph):
        spec, labels = self.spec, node_labels(graph)
        if any(label not in self.objects and label not in self.edge_nodes for label in labels.values()):
            return "unknown-node-sort"
        relation = defaultdict(Counter)
        for u, v, bond in edge_rows(graph):
            if bond in (spec["sourceBond"], spec["targetBond"]):
                if (labels[u] in self.edge_nodes) == (labels[v] in self.edge_nodes):
                    return "incidence-sort-conflict"
            elif bond in self.relations:
                if labels[u] in self.edge_nodes or labels[v] in self.edge_nodes:
                    return f"relation-edge-node-conflict:{bond}"
                pair = tuple(sorted((labels[u], labels[v])))
                if pair not in self.relations[bond]:
                    return f"forbidden-relation:{bond}:{'|'.join(pair)}"
                relation[u][bond, labels[v]] += 1
                relation[v][bond, labels[u]] += 1
            else:
                return f"unsupported-incidence:{bond}"
        ports = self._ports(graph)
        outgoing, incoming = Counter(), Counter()
        for edge, (sources, targets) in ports.items():
            if len(sources) > 1 or len(targets) > 1:
                return "nonfunctional-edge-port"
            outgoing.update(sources)
            if labels[edge] in spec["atMostOneIncomingTypes"]:
                incoming.update((target, labels[edge]) for target in targets)
            ends = spec.get("edgeEndpointTypes", {}).get(labels[edge])
            if ends is not None:
                if sources and labels[sources[0]] not in ends["source"]:
                    return f"forbidden-source-type:{labels[edge]}"
                if targets and labels[targets[0]] not in ends["target"]:
                    return f"forbidden-target-type:{labels[edge]}"
        for bound in spec.get("relationIncidentBounds", ()):
            for u, label in labels.items():
                if label == bound["center"] and sum(relation[u][bound["bond"], other]
                        for other in bound["neighbourLabels"]) > bound["maximum"]:
                    return f"relation-bound:{bound['name']}"
        if spec["atMostOneOutgoing"] and any(n > 1 for n in outgoing.values()):
            return "multiple-outgoing-edges"
        if any(n > 1 for n in incoming.values()):
            return "duplicate-incoming-edge-type"
        for edge, (_, targets) in ports.items():
            if labels[edge] in spec["forbidSuccessorAfterTypes"] and targets and outgoing[targets[0]]:
                return f"successor-after-{labels[edge]}"
        return self.directed_violation(graph)


def typed_incidence_constraints(specification: Mapping, vertex_labels, edge_labels):
    """MolFoundry incidence constraints for the local part of ``violation``, on
    graphs labelled from the given finite label sets."""
    from molfoundry import AnyOf, IncidenceConstraint, IncidenceCount
    spec = specification
    vertices, edges = sorted(vertex_labels), sorted(edge_labels)
    objects = [label for label in vertices if label in set(spec["objectLabels"])]
    edge_nodes = [label for label in vertices if label in set(spec["edgeNodeLabels"])]
    source, target = spec["sourceBond"], spec["targetBond"]
    relations = {bond: {tuple(sorted(pair)) for pair in pairs}
                 for bond, pairs in spec.get("relationBonds", {}).items() if bond in edges}
    constraints = []
    unknown = [label for label in vertices if label not in objects and label not in edge_nodes]
    if unknown:
        constraints.append(IncidenceConstraint(unknown, condition=AnyOf()))
    unsupported = [bond for bond in edges if bond not in (source, target) and bond not in relations]
    if unsupported:
        constraints.append(IncidenceConstraint(None, "<=", 0, edgeLabels=unsupported))
    if edge_nodes:
        if relations:
            constraints.append(IncidenceConstraint(edge_nodes, "<=", 0, edgeLabels=list(relations)))
        constraints += [
            IncidenceConstraint(edge_nodes, "<=", 0, edgeLabels=[source, target], nodeLabels=edge_nodes),
            IncidenceConstraint(edge_nodes, "<=", 1, edgeLabels=source),
            IncidenceConstraint(edge_nodes, "<=", 1, edgeLabels=target),
        ]
    if objects:
        constraints.append(IncidenceConstraint(objects, "<=", 0, edgeLabels=[source, target],
                                               nodeLabels=objects))
        for center in objects:
            for bond, pairs in relations.items():
                others = [o for o in vertices if tuple(sorted((center, o))) not in pairs]
                if others:
                    constraints.append(IncidenceConstraint(center, "<=", 0, nodeLabels=others,
                                                           edgeLabels=bond))
        if spec["atMostOneOutgoing"]:
            constraints.append(IncidenceConstraint(objects, "<=", 1, edgeLabels=source))
        for kind in spec["atMostOneIncomingTypes"]:
            if kind in edge_nodes:
                constraints.append(IncidenceConstraint(objects, "<=", 1, edgeLabels=target,
                                                       nodeLabels=kind))
        terminal = [kind for kind in spec["forbidSuccessorAfterTypes"] if kind in edge_nodes]
        if terminal:
            # An incoming terminal edge and an outgoing edge meet at the object.
            constraints.append(IncidenceConstraint(objects, condition=AnyOf(
                IncidenceCount("<=", 0, edgeLabels=target, nodeLabels=terminal),
                IncidenceCount("<=", 0, edgeLabels=source))))
    for kind, ends in spec.get("edgeEndpointTypes", {}).items():
        if kind in edge_nodes:
            for bond, role in ((source, "source"), (target, "target")):
                others = [o for o in objects if o not in ends[role]]
                if others:
                    constraints.append(IncidenceConstraint(kind, "<=", 0, edgeLabels=bond,
                                                           nodeLabels=others))
    for bound in spec.get("relationIncidentBounds", ()):
        if bound["center"] in objects and bound["bond"] in relations:
            constraints.append(IncidenceConstraint(bound["center"], "<=", bound["maximum"],
                edgeLabels=bound["bond"], nodeLabels=bound["neighbourLabels"]))
    return constraints


def grouped_typed_incidence_compositions(
    first,
    second,
    *,
    specification: Mapping,
    connected: bool = False,
    maximum: bool = False,
    include_empty: bool = True,
) -> tuple[IncidenceCompositionDiagnostics, Counter, dict[str, object]]:
    """Compose and reject violations of an incidence structural theory."""

    required = {
        "name", "objectLabels", "edgeNodeLabels", "sourceBond", "targetBond",
        "atMostOneOutgoing", "atMostOneIncomingTypes",
        "forbidSuccessorAfterTypes", "simpleDirected", "simpleUnderlying",
        "acyclicDirected",
    }
    optional = {"edgeEndpointTypes", "relationBonds", "relationIncidentBounds"}
    if not required <= set(specification) or set(specification) - required - optional:
        raise ValueError("typed-incidence specification has unexpected fields")
    object_labels = specification["objectLabels"]
    edge_labels = specification["edgeNodeLabels"]
    if (
        not object_labels
        or not edge_labels
        or set(object_labels) & set(edge_labels)
        or not set(specification["atMostOneIncomingTypes"]) <= set(edge_labels)
        or not set(specification["forbidSuccessorAfterTypes"]) <= set(edge_labels)
        or specification["sourceBond"] == specification["targetBond"]
        or type(specification["acyclicDirected"]) is not bool
    ):
        raise ValueError("invalid typed-incidence signature or local bounds")
    endpoint_types = specification.get("edgeEndpointTypes", {})
    if set(endpoint_types) - set(edge_labels):
        raise ValueError("endpoint typing refers to an unknown edge-node label")
    for edge, endpoints in endpoint_types.items():
        if set(endpoints) != {"source", "target"}:
            raise ValueError(f"invalid endpoint typing for {edge!r}")
        if (
            not endpoints["source"]
            or not endpoints["target"]
            or not set(endpoints["source"]) <= set(object_labels)
            or not set(endpoints["target"]) <= set(object_labels)
        ):
            raise ValueError(f"invalid endpoint labels for {edge!r}")
    relation_bonds = specification.get("relationBonds", {})
    if any(bond in {specification["sourceBond"], specification["targetBond"]}
           for bond in relation_bonds):
        raise ValueError("relation bonds must differ from incidence bonds")
    for bond, pairs in relation_bonds.items():
        if not isinstance(bond, str) or not bond or not pairs:
            raise ValueError("relation bonds need a nonempty label and pair set")
        if any(
            not isinstance(pair, (list, tuple))
            or len(pair) != 2
            or not set(pair) <= set(object_labels)
            for pair in pairs
        ):
            raise ValueError(f"invalid allowed relation pair for {bond!r}")
    for bound in specification.get("relationIncidentBounds", []):
        if set(bound) != {"center", "bond", "neighbourLabels", "maximum", "name"}:
            raise ValueError("invalid relation incident bound fields")
        if (
            bound["center"] not in object_labels
            or bound["bond"] not in relation_bonds
            or not bound["neighbourLabels"]
            or not set(bound["neighbourLabels"]) <= set(object_labels)
            or type(bound["maximum"]) is not int
            or bound["maximum"] < 0
            or not isinstance(bound["name"], str)
            or not bound["name"]
        ):
            raise ValueError("invalid relation incident bound")
    predicate = TypedIncidencePredicate(specification)
    global_part = specification["simpleDirected"] or specification["acyclicDirected"]
    return group_compositions(
        first,
        second,
        constraints=typed_incidence_constraints(specification, *rule_labels(first, second)),
        violation=predicate.directed_violation if global_part else None,
        connected=connected,
        maximum=maximum,
        include_empty=include_empty,
    )


@dataclass(frozen=True)
class IncidenceCommutator:
    coefficients: Mapping[str, int]
    representatives: Mapping[str, object]
    left: IncidenceCompositionDiagnostics
    right: IncidenceCompositionDiagnostics


def typed_incidence_commutator(
    observable,
    reaction,
    *,
    specification: Mapping,
    connected: bool = False,
    include_empty: bool = True,
) -> IncidenceCommutator:
    """Return the raw composite-rule expression ``reaction*O - O*reaction``."""

    left, left_counts, left_representatives = grouped_typed_incidence_compositions(
        observable,
        reaction,
        specification=specification,
        connected=connected,
        include_empty=include_empty,
    )
    right, right_counts, right_representatives = grouped_typed_incidence_compositions(
        reaction,
        observable,
        specification=specification,
        connected=connected,
        include_empty=include_empty,
    )
    coefficients = {
        key: right_counts[key] - left_counts[key]
        for key in set(left_counts) | set(right_counts)
        if right_counts[key] != left_counts[key]
    }
    return IncidenceCommutator(
        coefficients,
        left_representatives | right_representatives,
        left,
        right,
    )


def pair_commutator(specification: Mapping, *, include_empty: bool):
    """``typed_incidence_commutator(observable, reaction)`` under ``specification``,
    with unrestricted (not necessarily connected) overlaps: the ``pair_commutator``
    of the composition jump compiler and of observable discovery."""

    return partial(
        typed_incidence_commutator,
        specification=specification,
        connected=False,
        include_empty=include_empty,
    )


def dangling_certificate(specification: Mapping):
    """Certify DPO dangling of deleted vertices from the declared structural theory.

    A deleted vertex has no host edge outside the match if every incidence the
    theory admits at its sort is bounded and already saturated in ``I``: an
    edge-node has at most one source and one target port and no relation
    bonds; an object's outgoing ports, bounded incoming edge types and bounded
    relation bonds must reach their maxima.  Valid on every state satisfying
    ``specification``; anything else raises ``DanglingNotCertified``.
    """
    from mfobs.core.identity import DanglingNotCertified

    objects = set(specification["objectLabels"])
    edge_nodes = set(specification["edgeNodeLabels"])
    source_bond, target_bond = specification["sourceBond"], specification["targetBond"]
    endpoints = specification.get("edgeEndpointTypes", {})
    bounded_incoming = set(specification.get("atMostOneIncomingTypes", ()))
    relations = {
        bond: {label for pair in pairs for label in pair}
        for bond, pairs in specification.get("relationBonds", {}).items()
    }
    bounds = specification.get("relationIncidentBounds", ())

    def admits(edge_type: str, role: str, label: str) -> bool:
        allowed = endpoints.get(edge_type)
        return allowed is None or label in allowed[role]

    def certify(graph, deleted) -> None:
        labels = node_labels(graph)
        for vertex in deleted:
            label = labels[vertex]
            incident = [(v if u == vertex else u, bond)
                        for u, v, bond in edge_rows(graph) if vertex in (u, v)]
            bonds = Counter(bond for _, bond in incident)
            if label in edge_nodes:
                if bonds[source_bond] and bonds[target_bond]:
                    continue
                raise DanglingNotCertified(
                    f"deleted edge-node {vertex} ({label}) lacks a source or target port in I")
            if label not in objects:
                raise DanglingNotCertified(f"deleted vertex {vertex} has an undeclared sort {label}")
            ports = [(labels[other], bond) for other, bond in incident
                     if labels[other] in edge_nodes]
            if any(admits(t, "source", label) for t in edge_nodes):
                if not (specification.get("atMostOneOutgoing")
                        and any(bond == source_bond for _, bond in ports)):
                    raise DanglingNotCertified(
                        f"deleted object {vertex} ({label}) may have an unmatched outgoing edge")
            for edge_type in sorted(edge_nodes):
                if admits(edge_type, "target", label) and not (
                        edge_type in bounded_incoming
                        and (edge_type, target_bond) in ports):
                    raise DanglingNotCertified(
                        f"deleted object {vertex} ({label}) may have an unmatched incoming {edge_type}")
            for bond, sorts in relations.items():
                if label not in sorts:
                    continue
                neighbours = {m for pair in specification["relationBonds"][bond]
                              if label in pair for m in pair if m != label or pair.count(label) == 2}
                saturated = any(
                    bound["center"] == label and bound["bond"] == bond
                    and neighbours <= set(bound["neighbourLabels"])
                    and sum(1 for other, b in incident if b == bond
                            and labels[other] in bound["neighbourLabels"]) >= bound["maximum"]
                    for bound in bounds
                )
                if not saturated:
                    raise DanglingNotCertified(
                        f"deleted object {vertex} ({label}) may have an unmatched {bond} relation")

    return certify
