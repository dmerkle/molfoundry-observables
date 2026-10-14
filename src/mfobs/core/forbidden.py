"""Finite local forbidden-pattern specifications for restricted composition.

The TCS 2021 restricted-rewriting construction defines the pattern category by
excluding a finite family ``N`` of forbidden subgraphs.  For performance, this
prototype compiles a finite family of *local star patterns* into MolFoundry
incidence constraints (``incidence_constraints``), which the engine checks on
every overlap graph and uses to prune the overlap search.
``LocalPatternPredicate`` is the same check in Python.  The existing semantic
reducer remains the final authority for model-specific reachability identities.

A specification contains only monotone negative information: once one of the
listed local patterns is present in a graph, embedding that graph into a larger
chemical state cannot repair it. This applies to monic host extension, NOT
to extending a partially chosen overlap with further identifications.  Global cardinality facts (fixed catalyst populations, etc.) and positive completion constraints are deliberately not included here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations_with_replacement
from typing import Mapping


Pair = tuple[str, str]
Port = tuple[str, str]  # (edge label, neighbour atom label)
MAX_SAFE_INTEGER = 2**53 - 1


def _label(value, description):
    if not isinstance(value, str) or not value or '|' in value:
        raise ValueError(f'{description} must be a nonempty label without |')


def _nonnegative_integer(value, description):
    # The backend uses JavaScript Numbers. In particular bool is not an integer
    # in the input contract, and accepting a rounded float can change pruning.
    if type(value) is not int or not 0 <= value <= MAX_SAFE_INTEGER:
        raise ValueError(f'{description} must be a nonnegative safe integer')


def _pair(a: str, b: str) -> Pair:
    return tuple(sorted((a, b)))  # type: ignore[return-value]


@dataclass(frozen=True)
class IncidentBound:
    """Forbid more than ``maximum`` matching incident edges at one atom."""

    center: str
    edge: str
    neighbour: str
    maximum: int
    name: str

    def __post_init__(self):
        for field in ('center', 'edge', 'neighbour', 'name'):
            _label(getattr(self, field), f'incident bound {field}')
        _nonnegative_integer(self.maximum, 'incident bound maximum')


@dataclass(frozen=True)
class ForbiddenPatternSpecification:
    """Compiled finite connected forbidden-pattern family.

    ``allowed_bond_pairs`` is a compact complement representation of all
    forbidden two-atom/bond patterns over the finite signature.  ``max_valence``
    compiles the finite set of minimal over-valence star patterns.  The explicit
    incident bounds compile additional finite stars not implied by valence.
    """

    name: str
    atom_labels: tuple[str, ...]
    bond_orders: Mapping[str, int]
    allowed_bond_pairs: Mapping[str, tuple[Pair, ...]]
    max_valence: Mapping[str, int]
    incident_bounds: tuple[IncidentBound, ...] = ()
    strict_signature: bool = True

    def __post_init__(self):
        self.validate()
        # Both the executable checker and the explicit family use unordered
        # atom pairs. Canonicalize once so their certificates cannot disagree.
        object.__setattr__(self, 'allowed_bond_pairs', {
            edge: tuple(sorted({_pair(*pair) for pair in pairs}))
            for edge, pairs in self.allowed_bond_pairs.items()
        })

    def validate(self):
        _label(self.name, 'forbidden theory name')
        if not isinstance(self.atom_labels, (tuple, list)) or not self.atom_labels:
            raise ValueError('atomLabels must be a nonempty sequence')
        for label in self.atom_labels:
            _label(label, 'atom label')
        labels = set(self.atom_labels)
        if len(labels) != len(self.atom_labels):
            raise ValueError('atomLabels must be unique')
        if type(self.strict_signature) is not bool:
            raise ValueError('strictSignature must be boolean')
        for name, mapping in [('bondOrders', self.bond_orders),
                              ('allowedBondPairs', self.allowed_bond_pairs),
                              ('maxValence', self.max_valence)]:
            if not isinstance(mapping, Mapping):
                raise ValueError(f'{name} must be a mapping')
        if set(self.bond_orders) != set(self.allowed_bond_pairs):
            raise ValueError('bondOrders and allowedBondPairs must declare the same edges')
        for edge, order in self.bond_orders.items():
            _label(edge, 'bond label')
            _nonnegative_integer(order, f'bond order {edge}')
            pairs = self.allowed_bond_pairs[edge]
            if not isinstance(pairs, (tuple, list)):
                raise ValueError(f'allowedBondPairs[{edge}] must be a sequence')
            for pair in pairs:
                if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                    raise ValueError('allowed bond pairs must have two atom labels')
                if any(label not in labels for label in pair):
                    raise ValueError('allowed bond pair uses an undeclared atom label')
        if set(self.max_valence) - labels:
            raise ValueError('maxValence uses an undeclared atom label')
        for atom, maximum in self.max_valence.items():
            _nonnegative_integer(maximum, f'maximum valence of {atom}')
        if not isinstance(self.incident_bounds, (tuple, list)):
            raise ValueError('incidentBounds must be a sequence')
        for bound in self.incident_bounds:
            if not isinstance(bound, IncidentBound):
                raise ValueError('incidentBounds entries must be IncidentBound values')
            bound.__post_init__()
            if bound.center not in labels or bound.neighbour not in labels:
                raise ValueError('incident bound uses an undeclared atom label')
            if bound.edge not in self.bond_orders:
                raise ValueError('incident bound uses an undeclared edge label')

    def engine_payload(self) -> dict:
        self.validate()
        return {
            "name": self.name,
            "atomLabels": list(self.atom_labels),
            "bondOrders": dict(self.bond_orders),
            "allowedBondPairs": {
                label: [list(pair) for pair in pairs]
                for label, pairs in self.allowed_bond_pairs.items()
            },
            "maxValence": dict(self.max_valence),
            "incidentBounds": [asdict(bound) for bound in self.incident_bounds],
            "strictSignature": self.strict_signature,
        }

    def invalid_bond_patterns(self) -> list[dict]:
        """Expand the complement encoding to concrete two-node patterns."""

        all_pairs = list(combinations_with_replacement(self.atom_labels, 2))
        patterns: list[dict] = []
        for edge in self.bond_orders:
            allowed = {tuple(pair) for pair in self.allowed_bond_pairs.get(edge, ())}
            for a, b in all_pairs:
                pair = _pair(a, b)
                if pair not in allowed:
                    patterns.append(
                        {
                            "kind": "invalid-bond-signature",
                            "nodes": [a, b],
                            "edge": edge,
                        }
                    )
        return patterns

    def minimal_overvalence_patterns(self) -> list[dict]:
        """Enumerate minimal forbidden local valence stars.

        Zero-order relations (if a schema has one) do not participate
        in covalent valence and are handled by explicit incident bounds.
        """

        incident_by_center: dict[str, list[Port]] = {
            label: [] for label in self.atom_labels
        }
        for edge, pairs in self.allowed_bond_pairs.items():
            order = self.bond_orders[edge]
            if order <= 0:
                continue
            for a, b in pairs:
                incident_by_center[a].append((edge, b))
                if a != b:
                    incident_by_center[b].append((edge, a))
                else:
                    incident_by_center[b].append((edge, a))

        patterns: list[dict] = []
        for center, maximum in self.max_valence.items():
            ports = sorted(set(incident_by_center.get(center, [])))
            # At least one positive-order edge is added each time, so maximum+1
            # incident edges suffice to find all minimal excess stars.
            for size in range(1, maximum + 2):
                for combo in combinations_with_replacement(ports, size):
                    orders = [self.bond_orders[edge] for edge, _ in combo]
                    total = sum(orders)
                    if total <= maximum:
                        continue
                    if any(total - order > maximum for order in orders):
                        continue
                    patterns.append(
                        {
                            "kind": "over-valence-star",
                            "center": center,
                            "ports": [
                                {"edge": edge, "neighbour": neighbour}
                                for edge, neighbour in combo
                            ],
                            "valence": total,
                            "maximum": maximum,
                        }
                    )
        return patterns

    def certificate(self) -> dict:
        self.validate()
        bond_patterns = self.invalid_bond_patterns()
        valence_patterns = self.minimal_overvalence_patterns()
        incidence_patterns = [
            {
                "kind": "incident-bound-star",
                "center": bound.center,
                "edge": bound.edge,
                "neighbour": bound.neighbour,
                "required_distinct_neighbours": bound.maximum + 1,
                "name": bound.name,
            }
            for bound in self.incident_bounds
        ]
        return {
            "name": self.name,
            "semantics": (
                "finite connected negative-pattern family compiled into a local "
                "overlap-pushout checker"
            ),
            "scope": (
                "negative local pattern constraints only; excludes global "
                "cardinality and positive state-completion invariants"
            ),
            "concrete_invalid_bond_patterns": bond_patterns,
            "minimal_overvalence_star_patterns": valence_patterns,
            "incident_bound_star_patterns": incidence_patterns,
            "counts": {
                "invalid_bond_patterns": len(bond_patterns),
                "minimal_overvalence_patterns": len(valence_patterns),
                "incident_bound_patterns": len(incidence_patterns),
                "total_compiled_patterns": (
                    len(bond_patterns)
                    + len(valence_patterns)
                    + len(incidence_patterns)
                ),
            },
        }


class LocalPatternPredicate:
    """Reject graphs containing a forbidden bond, star or neighbourhood."""
    def __init__(self, spec):
        self.spec = spec
        self.known = set(spec['atomLabels'])
        self.allowed = {bond: {tuple(sorted(pair)) for pair in pairs}
                        for bond, pairs in spec['allowedBondPairs'].items()}

    def violation(self, graph):
        from collections import Counter, defaultdict
        from mfobs.core.graph import node_labels, edge_rows
        spec, labels = self.spec, node_labels(graph)
        if spec['strictSignature'] and set(labels.values()) - self.known:
            return 'unknown-atom-label'
        valence, incident = Counter(), defaultdict(Counter)
        for u, v, bond in edge_rows(graph):
            if bond not in self.allowed:
                return f'unsupported-bond:{bond}'
            pair = tuple(sorted((labels[u], labels[v])))
            if pair not in self.allowed[bond]:
                return f'forbidden-bond:{bond}:{"|".join(pair)}'
            valence[u] += spec['bondOrders'][bond]
            valence[v] += spec['bondOrders'][bond]
            incident[u][bond, labels[v]] += 1
            incident[v][bond, labels[u]] += 1
        for u, label in labels.items():
            maximum = spec['maxValence'].get(label)
            if maximum is not None and valence[u] > maximum:
                return f'forbidden-valence:{label}:{valence[u]}>{maximum}'
        for bound in spec['incidentBounds']:
            for u, label in labels.items():
                if label == bound['center'] and incident[u][bound['edge'], bound['neighbour']] > bound['maximum']:
                    return f'forbidden-star:{bound["name"]}'
        profiles = spec.get('allowedIncidentProfiles')
        if profiles is not None:
            for u, label in labels.items():
                if not any(all(n <= profile.get(bond+'|'+other, 0)
                               for (bond, other), n in incident[u].items())
                           for profile in profiles.get(label, [])):
                    return f'forbidden-neighbourhood:{label}'
        return None


def incidence_constraints(spec, vertex_labels, edge_labels):
    """MolFoundry incidence constraints rejecting exactly the graphs that
    ``LocalPatternPredicate(spec).violation`` rejects, among graphs labelled
    from the given finite label sets (see ``composition.rule_labels``).

    Every check of the predicate is a condition on the edges at one vertex. An
    unknown label or a label without an allowed profile is excluded outright
    (an empty ``AnyOf`` never holds); complements such as "any other neighbour"
    are taken within the label sets.
    """
    from molfoundry import AllOf, AnyOf, IncidenceConstraint, IncidenceCount
    vertices, edges = sorted(vertex_labels), sorted(edge_labels)
    allowed = {bond: {tuple(sorted(pair)) for pair in pairs}
               for bond, pairs in spec['allowedBondPairs'].items()}
    constraints = []
    unknown = [label for label in vertices if label not in set(spec['atomLabels'])]
    if spec['strictSignature'] and unknown:
        constraints.append(IncidenceConstraint(unknown, condition=AnyOf()))
    unsupported = [bond for bond in edges if bond not in allowed]
    if unsupported:
        constraints.append(IncidenceConstraint(None, '<=', 0, edgeLabels=unsupported))
    for center in vertices:
        for bond in edges:
            if bond in allowed:
                others = [o for o in vertices if tuple(sorted((center, o))) not in allowed[bond]]
                if others:
                    constraints.append(IncidenceConstraint(center, '<=', 0, nodeLabels=others,
                                                           edgeLabels=bond))
        if center in spec['maxValence']:
            orders = {bond: spec['bondOrders'][bond] for bond in edges if bond in allowed}
            constraints.append(IncidenceConstraint(center, '<=', spec['maxValence'][center],
                                                   edgeWeights=orders))
    for bound in spec['incidentBounds']:
        if bound['center'] in vertices:
            constraints.append(IncidenceConstraint(bound['center'], '<=', bound['maximum'],
                nodeLabels=bound['neighbour'], edgeLabels=bound['edge']))
    profiles = spec.get('allowedIncidentProfiles')
    if profiles is not None:
        for center in vertices:
            options = []
            for profile in profiles.get(center, []):
                ports = {tuple(key.split('|', 1)): n for key, n in profile.items()}
                counts = [IncidenceCount('<=', n, nodeLabels=other, edgeLabels=bond)
                          for (bond, other), n in sorted(ports.items())
                          if bond in edges and other in vertices]
                for bond in edges:
                    others = [o for o in vertices if (bond, o) not in ports]
                    if others:
                        counts.append(IncidenceCount('<=', 0, nodeLabels=others, edgeLabels=bond))
                options.append(AllOf(*counts))
            constraints.append(IncidenceConstraint(center, condition=AnyOf(*options)))
    return constraints
