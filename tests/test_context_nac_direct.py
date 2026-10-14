"""Concrete differential audit, independent of the Boolean reducer.

Count exported positive patterns with NetworkX on actual DG products. Compare
complete jump rates with generated events and a separate four-site reference.
Check transported conditions term by term before commutator cancellation.
Finite samples test the implementation; all-tail scope rests on locality proofs.
"""
import mfobs
from collections import Counter, defaultdict

import molfoundry as mf
import networkx as nx
import pytest
import sympy as sp

from mfobs import Chlorination
from mfobs.core.engine import commutator, dpo_jump_closure
from mfobs.core.gml import graph_gml
from mfobs.models.chlorination.inputs import methane
from mfobs.checks.chlorination import _molecule, _hcl


def counts(molecule, patterns):
    """Literal positive embeddings, without the synthesis reducer or its matches."""
    host = nx.Graph()
    host.add_nodes_from((v.id, {"label": v.stringLabel}) for v in molecule.vertices)
    host.add_edges_from((e.source.id, e.target.id, {"label": e.stringLabel})
                        for e in molecule.edges)
    result = {}
    for name, pattern in patterns.items():
        query = nx.Graph()
        query.add_nodes_from((v.id, {"label": v.leftLabel})
                             for v in pattern.rule.vertices if v.leftLabel is not None)
        query.add_edges_from((e.source.id, e.target.id, {"label": e.leftLabel})
                             for e in pattern.rule.edges if e.leftLabel is not None)
        matcher = nx.isomorphism.GraphMatcher(
            host, query, node_match=lambda a, b: a == b, edge_match=lambda a, b: a == b)
        result[name] = sum(1 for _ in matcher.subgraph_monomorphisms_iter())
    return result


def applied(components, rule):
    dg = mf.DG(graphDatabase=components)
    with dg.build() as builder:
        edges = builder.apply(components, rule, mode=mf.ApplyMode.UseAny, verbosity=0)
    return dg, edges


@pytest.fixture(scope="module", params=["guarded", "ablated"])
def analyzed(request):
    backend = Chlorination(request.param)
    return request.param, backend, mfobs.analyze(backend.specification(), backend=backend).record()


def test_joint_law_and_transported_terms_against_public_dg(analyzed):
    variant, backend, result = analyzed
    patterns = backend.patterns
    names = tuple(patterns)
    # Check each literal projected term, not only its summed generator image.
    activity = dpo_jump_closure(backend.rule)
    terms = {activity.canonKey: activity}
    for pattern in patterns.values():
        composed = commutator(pattern.rule, backend.rule)
        terms.update(composed.observable_representatives)
    expressions = {key: backend.reducer.reduce(term).expression for key, term in terms.items()}
    joint_checks = transport_checks = 0
    for n in (0, 1, 7):
        values = {mask: counts(_molecule(mask, n), patterns) for mask in range(16)}
        for mask in range(16):
            molecule, reagent = _molecule(mask, n), _hcl()
            q = 4 - mask.bit_count()
            components = [molecule] + ([reagent] if q else [])
            initial = values[mask]
            substitution = {sp.Symbol(name): value for name, value in initial.items()}
            dg, edges = applied(components, backend.rule)
            actual = defaultdict(int)
            for edge in edges:
                # DG counts include graph symmetries. The distinct reactants
                # consume one head and one HCl, giving inventory factor q.
                assert Counter(v.graph.canonKey for v in edge.sources) == Counter(
                    [molecule.canonKey, reagent.canonKey])
                after = dict.fromkeys(names, 0)
                for vertex in edge.targets:
                    for name, value in counts(vertex.graph, patterns).items():
                        after[name] += value
                delta = tuple(after[name] - initial[name] for name in names)
                actual[delta] += edge.numMatches(backend.rule) * q
            expected = defaultdict(int)
            for event in result["events"]:
                rate = sp.sympify(event["activity"]["text"]).subs(substitution)
                assert rate.is_Integer and rate >= 0
                if rate:
                    expected[tuple(event["delta"].get(name, 0) for name in names)] += int(rate)
            # Separate concrete path-four reference; no synthesis enablement helper.
            reference = defaultdict(int)
            for site in range(4):
                if mask & (1 << site):
                    continue
                if variant == "guarded" and any(
                        mask & (1 << neighbour) for neighbour in (site - 1, site + 1)
                        if 0 <= neighbour < 4):
                    continue
                after = values[mask | (1 << site)]
                if q:
                    reference[tuple(after[name] - initial[name] for name in names)] += q
            assert actual == expected == reference, (variant, n, mask)
            joint_checks += 1
            for key, term in terms.items():
                dg, edges = applied(components, term)
                total = 0
                for edge in edges:
                    sources = Counter(v.graph.canonKey for v in edge.sources)
                    assert all(number == 1 for number in sources.values())
                    assert set(sources) <= {molecule.canonKey, reagent.canonKey}
                    total += edge.numMatches(term) * (q if reagent.canonKey in sources else 1)
                assert total == expressions[key].subs(substitution), (
                    variant, n, mask, term.getGMLString())
                transport_checks += 1
    assert joint_checks == 48
    assert transport_checks == 48 * len(terms)
    if variant == "guarded":
        assert any(term.matchCondition is not None for term in terms.values())
    else:
        assert all(term.matchCondition is None for term in terms.values())


def test_growth_rate_and_zero_observable_jump(analyzed):
    _, backend, _ = analyzed
    reagent = mf.Graph.fromGMLString(graph_gml(methane()), add=False)
    for n in (0, 1, 7):
        for mask in range(16):
            molecule = _molecule(mask, n)
            before = counts(molecule, backend.patterns)
            for inventory in (0, 3):
                components = [molecule] + ([reagent] if inventory else [])
                dg, edges = applied(components, backend.growth_rule)
                rate = 0
                for edge in edges:
                    assert Counter(v.graph.canonKey for v in edge.sources) == Counter(
                        [molecule.canonKey, reagent.canonKey])
                    after = dict.fromkeys(before, 0)
                    for vertex in edge.targets:
                        for name, value in counts(vertex.graph, backend.patterns).items():
                            after[name] += value
                    assert after == before
                    rate += edge.numMatches(backend.growth_rule) * backend.growth_rate_factor * inventory
                assert rate == inventory  # k_grow = 1, including the 1/144 factor.


def test_ablation_removes_only_the_generic_condition():
    guarded = Chlorination()
    ablated = Chlorination("ablated")
    assert guarded.rule.matchCondition is not None
    assert ablated.rule.matchCondition is None
    assert guarded.rule.getGMLString().split("\tcondition [")[0] + "]" == ablated.rule.getGMLString()
    assert guarded.growth_rule.getGMLString() == ablated.growth_rule.getGMLString()
    assert guarded.chlorine_total == ablated.chlorine_total == 4
    assert Chlorination().query_gmls() == Chlorination("ablated").query_gmls()


def test_nac_analysis_enumerates_no_matches(monkeypatch):
    from mfobs.core import molfoundry_execution

    def forbidden(*args, **kwargs):
        raise AssertionError("NAC analysis enumerated matches through the execution adapter")

    monkeypatch.setattr(molfoundry_execution, "matches", forbidden)
    backend1 = Chlorination()
    assert mfobs.analyze(backend1.specification(), backend=backend1).record()["status"] == "PASS"
