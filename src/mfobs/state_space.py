"""Inspectable state families and their shared local chemical restrictions."""
from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from importlib.resources import files
import json
from math import lcm
from dataclasses import dataclass
from types import MappingProxyType


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({k: _freeze(v) for k, v in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(v) for v in value)
    return value


def _record(value):
    if isinstance(value, Mapping):
        return {k: _record(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_record(v) for v in value]
    return value


@dataclass(frozen=True)
class StateSpace:
    """A backend's mathematical domain, not a membership decision procedure.

    Chemical validation checks the shared signature and upper bounds on rule
    sides or graph patterns. The backend additionally enforces its own family.
    """
    kind: str
    definition: Mapping
    chemistry: Mapping | None = None

    def __post_init__(self):
        if not isinstance(self.kind, str) or not self.kind:
            raise ValueError("state-space kind must be a nonempty string")
        if not isinstance(self.definition, Mapping):
            raise TypeError("state-space definition must be a mapping")
        object.__setattr__(self, "definition", _freeze(self.definition))
        object.__setattr__(self, "chemistry", _freeze(self.chemistry))

    def record(self):
        result = {"kind": self.kind, "definition": _record(self.definition)}
        if self.chemistry is not None:
            result["chemistry"] = _record(self.chemistry)
        return result

    def require(self, expected):
        """Refuse changes that a specialized backend has not implemented."""
        if self.record() != expected.record():
            raise ValueError(f"unsupported state space for {expected.kind}")
        return self

    def chemical_constraints(self):
        if self.chemistry is None:
            raise ValueError("this is not a chemical state space")
        return _record(self.chemistry["constraints"])

    def validate_graph(self, graph):
        """Check a partial graph's local chemistry, not full family membership."""
        if self.chemistry is None:
            return
        from mfobs.core.forbidden import LocalPatternPredicate
        reason = LocalPatternPredicate(self.chemical_constraints()).violation(graph)
        if reason:
            raise ValueError(f"{self.kind}: {reason}")

    def validate_rule(self, rule):
        """Check both actual MolFoundry rule sides without changing its span."""
        if self.chemistry is None:
            return
        from mfobs.core.rules import rule_graph
        self.validate_graph(rule_graph(rule, 'left'))
        self.validate_graph(rule_graph(rule, 'right'))


@lru_cache(None)
def covalent_theory(coordination=False, feed_symbols=()):
    path = files("mfobs").joinpath("data/ordinary-covalent-v1.json")
    data = json.loads(path.read_text(encoding="utf-8"))
    return chemical_theory(data, "mfobs/data/ordinary-covalent-v1.json", coordination, feed_symbols)


def chemical_theory(data, source, coordination=False, feed_symbols=()):
    from mfobs.specifications import validate_chemical_specification
    validate_chemical_specification(data)
    orders = {k: Fraction(v["order"]) for k, v in data["bondTypes"].items()}
    bounds = {k: Fraction(v["maxValence"]) for k, v in data["atomTypes"].items()}
    pairs = {k: list(v) for k, v in data["allowedBondPairs"].items()}
    atoms = list(data["atomTypes"])
    extensions = {}
    if coordination:
        atoms += ["Pd", "P"]
        orders["coord"] = Fraction(0)
        pairs["coord"] = [["Pd", "P"], ["Pd", "O"], ["Pd", "N"]]
        pairs["-"] += [["P", "C"]]
        bounds["P"] = Fraction(3)
        extensions["coordination"] = {
            "atoms": ["Pd", "P"], "bond": "coord", "bond_order": 0,
            "pairs": pairs["coord"], "covalent_pairs": [["P", "C"]],
            "P_max_covalent_valence": 3,
            "Pd": "coordination only; no coordination-number bound",
        }
    if feed_symbols:
        markers = ["Feed-" + a for a in feed_symbols]
        atoms += markers
        extensions["boundary_markers"] = markers
        # No allowed bond pair contains a marker: each remains isolated.
    scale = lcm(*(v.denominator for v in (*orders.values(), *bounds.values())))
    return {
        "source": source,
        "scope": data["scope"], "extensions": extensions,
        "constraints": {
            "name": data["name"], "atomLabels": atoms,
            "bondOrders": {k: int(v * scale) for k, v in orders.items()},
            "allowedBondPairs": pairs,
            "maxValence": {k: int(v * scale) for k, v in bounds.items()},
            "incidentBounds": [], "strictSignature": True,
        },
    }


def alphabet_symbols(alphabet):
    alphabet = tuple(alphabet)
    if not alphabet or len(set(alphabet)) != len(alphabet) or any(
        not isinstance(a, str) or len(a) != 1 or not a.isascii() or not a.isalpha()
        for a in alphabet
    ):
        raise ValueError("alphabet must contain distinct single ASCII letters")
    return alphabet


def covalent_compositions():
    return StateSpace("local-covalent-compositions", {
        "graphs": "finite chemical graph patterns obeying the local signature and upper bounds",
        "scope": "local composition restrictions, not a full molecular state-family theorem",
    }, covalent_theory())


def from_manifest(inputs):
    definition = inputs.read("state_language")
    refs = inputs.data.get("specifications")
    if refs:
        data = json.loads(inputs.file(refs["chemistry"]).read_text(encoding="utf-8"))
        chemistry = chemical_theory(data, refs["chemistry"])
    else:
        chemistry = covalent_theory() if any(r["kind"] == "chemical" for r in inputs.data["rules"]) else None
    kind = definition.get("kind", definition.get("name", inputs.data["id"]))
    return StateSpace(kind, definition, chemistry)


def simple_graphs(*, vertex_labels, edge_labels):
    """All finite simple undirected graphs over the given label signature."""
    vertices, edges = tuple(sorted(set(vertex_labels))), tuple(sorted(set(edge_labels)))
    if not vertices or not edges:
        raise ValueError("simple graphs require nonempty vertex_labels and edge_labels")
    return StateSpace("finite-simple-graphs", {"vertex_labels": vertices, "edge_labels": edges})
