"""Sound descendant pruning must allow future overlaps to repair current ones."""
from mfobs.core.graph import make_graph
from mfobs.core.rules import make_rule
import random

from mfobs.core import engine
from mfobs.core.composition import group_compositions, rule_labels
from mfobs.core.forbidden import LocalPatternPredicate, incidence_constraints
from mfobs.core.typed_incidence import TypedIncidencePredicate, typed_incidence_constraints


from mfobs.core.gml import rule_gml


def test_incident_profiles_preserve_full_overlap_repair_and_multiplicity():
    graph = make_graph({0: "C", 1: "H"}, [(0, 1, "-")])
    bond = engine.mf.Rule.fromGMLString(
        rule_gml(make_rule("bond", graph, graph)), add=False)
    spec = {
        "atomLabels": ["C", "H"], "allowedBondPairs": {"-": [["C", "H"]]},
        "bondOrders": {"-": 1}, "maxValence": {"C": 1, "H": 1},
        "incidentBounds": [], "strictSignature": True,
    }
    expected = engine.grouped_restricted_common_compositions(
        bond, bond, forbidden_specification=spec)
    spec["allowedIncidentProfiles"] = {"C": [{"-|H": 1}], "H": [{"-|C": 1}]}
    actual = engine.grouped_restricted_common_compositions(
        bond, bond, forbidden_specification=spec)
    assert actual[1] == expected[1]
    assert sum(actual[1].values()) == 2  # full and empty; singleton can be repaired


def test_representatives_are_rebuilt_from_gml_in_their_class():
    graph = make_graph({0: "C", 1: "H", 2: "H"}, [(0, 1, "-"), (0, 2, "-")])
    rule = engine.mf.Rule.fromGMLString(
        rule_gml(make_rule("CH2", graph, graph)), add=False)
    _, left_counts, left = group_compositions(rule, rule)
    _, right_counts, right = group_compositions(rule, rule, include_empty=False)
    assert set(left) == set(left_counts) and set(right) == set(right_counts)
    merged = left | right
    assert set(merged) == set(left_counts) | set(right_counts)
    for key in merged:
        assert merged[key].canonKey == key
        assert merged[key].getGMLString() == merged.gml(key)


def _random_graphs(labels, bonds, count, seed, size=5):
    rng = random.Random(seed)
    for _ in range(count):
        n = rng.randint(1, size)
        nodes = {i: rng.choice(labels) for i in range(n)}
        yield make_graph(nodes, [(u, v, rng.choice(bonds)) for u in range(n)
                                 for v in range(u + 1, n) if rng.random() < 0.5])


def _accepts(graph, checks):
    """Whether ``graph`` survives as the only overlap graph of (empty -> graph) * empty."""
    make = engine.mf.Rule.fromGMLString(
        rule_gml(make_rule("G", make_graph({}), graph)), add=False)
    empty = engine.mf.Rule.fromGMLString('rule [ ruleID "empty" ]', add=False)
    return bool(group_compositions(make, empty, **checks(make, empty))[1])


def test_compiled_constraints_reject_exactly_the_local_pattern_violations():
    spec = {
        "atomLabels": ["C", "H", "O"], "strictSignature": True,
        "allowedBondPairs": {"-": [["C", "C"], ["C", "H"], ["C", "O"], ["H", "O"]], "=": [["C", "O"]]},
        "bondOrders": {"-": 1, "=": 2}, "maxValence": {"C": 3, "H": 1, "O": 2},
        "incidentBounds": [{"center": "C", "edge": "-", "neighbour": "H", "maximum": 1, "name": "one-h"}],
        "allowedIncidentProfiles": {"C": [{"-|H": 1, "-|C": 2}, {"=|O": 1, "-|C": 1}],
                                    "H": [{"-|C": 1}, {"-|O": 1}], "O": [{"=|C": 1}, {"-|H": 1, "-|C": 1}]},
    }
    loose = {**spec, "strictSignature": False}
    del loose["allowedIncidentProfiles"]
    for theory in (spec, loose):
        predicate = LocalPatternPredicate(theory)
        outcomes = set()
        for graph in _random_graphs(["C", "H", "O", "X"], ["-", "=", "#"], 150, seed=7):
            expected = predicate.violation(graph) is None
            outcomes.add(expected)
            assert _accepts(graph, lambda a, b: dict(
                constraints=incidence_constraints(theory, *rule_labels(a, b)))) == expected, graph
        assert outcomes == {True, False}


def test_typed_incidence_constraints_and_directed_predicate_match_the_reference():
    spec = {
        "name": "test", "objectLabels": ["P", "Q", "E"], "edgeNodeLabels": ["N", "T"],
        "sourceBond": "src", "targetBond": "tgt", "atMostOneOutgoing": True,
        "atMostOneIncomingTypes": ["N", "T"], "forbidSuccessorAfterTypes": ["T"],
        "simpleDirected": True, "simpleUnderlying": True, "acyclicDirected": True,
        "edgeEndpointTypes": {"N": {"source": ["P", "Q"], "target": ["P", "Q"]},
                              "T": {"source": ["P", "Q"], "target": ["E"]}},
        "relationBonds": {"rel": [["P", "Q"]]},
        "relationIncidentBounds": [{"center": "P", "bond": "rel", "neighbourLabels": ["Q"],
                                    "maximum": 1, "name": "one-rel"}],
    }
    rng = random.Random(11)

    def directed(count):
        # Encoded functional digraphs, mostly well typed, so that cycles and
        # (anti)parallel arcs are what decides.
        for _ in range(count):
            k = rng.randint(2, 4)
            nodes, edges = {i: rng.choice("PPQQE") for i in range(k)}, []
            targets = rng.sample(range(k), k)
            for source, target in enumerate(targets):
                if source != target and rng.random() < 0.85:
                    edge = len(nodes)
                    nodes[edge] = "N" if rng.random() < 0.85 else "T"
                    edges += [(edge, source, "src"), (edge, target, "tgt")]
            if rng.random() < 0.3:
                edges.append((0, 1, "rel"))
            yield make_graph(nodes, edges)

    predicate = TypedIncidencePredicate(spec)
    reasons = set()
    graphs = [*_random_graphs(["P", "Q", "E", "N", "T", "X"], ["src", "tgt", "rel", "bad"], 120, seed=3),
              *directed(120)]
    for graph in graphs:
        reason = predicate.violation(graph)
        reasons.add(reason)
        assert _accepts(graph, lambda a, b: dict(
            constraints=typed_incidence_constraints(spec, *rule_labels(a, b)),
            violation=predicate.directed_violation)) == (reason is None), (graph, reason)
    assert {None, "directed-cycle", "parallel-directed-edge"} <= reasons
