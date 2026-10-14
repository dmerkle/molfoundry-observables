#!/usr/bin/env python3
"""Independent analytical and full-graph checks for declarative gated paths.

This module is evaluation code, never imported by the synthesis adapter. It
does not use ChainLanguage, Reducer, stored reference equations, or a proposed
target basis. In particular, it checks whatever graph patterns were discovered.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import perm
from pathlib import Path
import sys
from time import perf_counter

import networkx as nx
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mfobs.core.conditions import transitions, engine_no_edges, engine_conditions
from mfobs.core.rules import rule_graph
from mfobs.core.engine import mf
from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows
from mfobs.inputs import load
from mfobs.checks.observable_expressions import parse_expressions, count

n, m, x = sp.symbols("n m x", integer=True)


def host(path_edges: int, tokens: int, opened: bool) -> mf.Graph:
    """Direct constructor, independent of the supplied block grammar."""
    nodes = {i: "A" for i in range(path_edges + 1)}
    edges = [(i, i + 1, "-") for i in range(path_edges)]
    u = path_edges + 1
    nodes.update({u: "U", u + 1: "W", u + 2: "V"})
    edges.extend([(u, u + 1, "-"), (u + 1, u + 2, "-")])
    if not opened:
        edges.append((u, u + 2, "-"))
    nodes.update({u + 3 + i: "T" for i in range(tokens)})
    return make_graph(nodes, edges)


def coordinates(graph: mf.Graph) -> tuple[int, int, int]:
    """Recognize the invariant from actual graph shape, not block provenance."""
    labels = node_labels(graph)
    path = {v for v, label in node_rows(graph) if label == "A"}
    control = {label: v for v, label in node_rows(graph) if label in {"U", "V", "W"}}
    assert len(control) == 3
    assert all(sum(label == want for _, label in node_rows(graph)) == 1 for want in control)
    assert all(label in {"A", "T", "U", "V", "W"} for label in labels.values())
    h = nx.Graph()
    h.add_nodes_from(path)
    h.add_edges_from((a, b) for a, b, _ in edge_rows(graph) if a in path and b in path)
    assert path and nx.is_tree(h) and max(dict(h.degree()).values()) <= 2
    expected_control = {(min(control[a], control[b]), max(control[a], control[b]))
                        for a, b in [("U", "W"), ("W", "V")]}
    actual_control = set()
    for a, b, label in edge_rows(graph):
        assert label == "-"
        if a in path or b in path:
            assert a in path and b in path
        else:
            assert labels[a] != "T" and labels[b] != "T"
            actual_control.add((a, b))
    optional = (min(control["U"], control["V"]), max(control["U"], control["V"]))
    assert actual_control in (expected_control, expected_control | {optional})
    return len(path) - 1, sum(label == "T" for label in labels.values()), int(optional not in actual_control)


def expression(form: dict) -> sp.Expr:
    result = 0
    for key, value in form.items():
        kind, degree, factorial_degree, binary_degree = key.split("|")
        assert kind == "power", "Unexpected finite-length exception in this path model"
        factorial = sp.prod(m - i for i in range(int(factorial_degree)))
        result += sp.Rational(value) * n**int(degree) * factorial * x**int(binary_degree)
    return sp.expand(result)


def boolean_reduce(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr.subs(x, 0) + x * (expr.subs(x, 1) - expr.subs(x, 0)))


def analytic_cyclic_rank(degree: int) -> dict:
    """Exact specialization diagnostic; NOT used to guide discovery.

    Rows are in ordinary powers n^j and n^j*x, with rates (2,3,5). A rank
    deficiency in the generated graph basis is reported, never hidden. For
    degree>=2 the invariant ideal n(n+1) explains the all-parameter upper bound.
    """
    assert degree >= 2
    matrix = sp.zeros(2 * degree)
    for j in range(1, degree + 1):
        for h in range(j):
            matrix[j - 1, degree + h] = matrix[degree + j - 1, degree + h] = 2 * sp.binomial(j, h)
        matrix[degree + j - 1, j - 1] = 5
        matrix[degree + j - 1, degree + j - 1] -= 8
    polynomial = sp.Poly(sp.prod(n + 1 - i for i in range(degree)), n)
    row = sp.Matrix([[polynomial.nth(i) for i in range(1, degree + 1)] + [0] * degree])
    rows = []
    for _ in range(2 * degree):
        rows.append(row)
        row = row * matrix
    rank = sp.Matrix.vstack(*rows).rank()
    return {"method": "exact-integer-Krylov-specialization-of-independent-generator",
            "degree": degree, "parameters": {"lambda_g": 2, "a": 3, "b": 5},
            "rank": rank, "polynomial_space_dimension": 2 * degree,
            "structural_upper_bound": 2 * degree - 2,
            "reason": "Seed belongs to invariant ideal n(n+1); rising-factorial degrees 2..k span a 2(k-1)-dimensional closure. Graph-representative discovery is not minimal."}


def oracle(f: sp.Expr, reaction: str, gated: bool, feedback: bool = False) -> sp.Expr:
    """Birth generator derived from oriented edge count 2n, not composition."""
    if reaction == "subdivide":
        return boolean_reduce(n * (x if gated else 1) * (f.subs(n, n + 1) - f))
    if reaction == "close":
        return boolean_reduce((n + 1 if feedback else 1) * x * (f.subs(x, 0) - f))
    if reaction == "open":
        return boolean_reduce((1 - x) * (f.subs(x, 1) - f))
    raise AssertionError(reaction)


def verify(inputs, result, *, max_n: int = 8, gated: bool = True,
           feedback: bool = False, partial: bool = False) -> dict:
    start = perf_counter()
    if partial:
        assert not result["closed"] and result["matrix"] is None and result["pending"]
        assert all(item["reason"] == "expansion budget exhausted" for item in result["unknown"])
    else:
        assert result["closed"] and not result["unknown"]
    assert inputs.definition_digest() == result["input_definition_sha256"]
    forms = {name: expression(form) for name, form in result["basis_counting_functions"].items()}
    observables = parse_expressions(result)
    rules = inputs.engine_rules()
    symbolic = 0
    for name, form in forms.items():
        if name not in result["equations_by_reaction"]:
            continue
        for reaction in rules:
            row = result["equations_by_reaction"][name][reaction]
            derived = sum(sp.Rational(c) * forms[target] for target, c in row.items())
            assert boolean_reduce(derived - oracle(form, reaction, gated, feedback)) == 0, (name, reaction)
            symbolic += 1
    count_checks = generator_checks = products = states = propensity_checks = 0
    for edges in range(max_n + 1):
        for tokens in (0, 2):
            for opened in (False, True):
                graph = host(edges, tokens, opened)
                assert coordinates(graph) == (edges, tokens, int(opened))
                states += 1
                before = {name: count(p, graph)
                          for name, p in observables.items()}
                for name, value in before.items():
                    assert value == forms[name].subs({n: edges, m: tokens, x: int(opened)}), (name, edges, opened)
                    count_checks += 1
                for declaration in inputs.data["rules"]:
                    reaction = declaration["id"]
                    targets = transitions(rules[reaction], graph)
                    expected_raw = (2 * edges * (int(opened) if gated else 1) if reaction == "subdivide"
                                    else (edges + 1 if feedback else 1) * int(opened) if reaction == "close"
                                    else int(not opened))
                    assert sum(targets.values()) == expected_raw
                    propensity_checks += 1
                    changes = {name: Fraction() for name in before}
                    for product, multiplicity in targets.items():
                        products += 1
                        expected = ((edges + 1, tokens, int(opened)) if reaction == "subdivide"
                                    else (edges, tokens, 0) if reaction == "close" else (edges, tokens, 1))
                        assert coordinates(product) == expected, (reaction, expected)
                        for name, pattern in observables.items():
                            changes[name] += (multiplicity * Fraction(declaration["per_match"]) *
                                (count(pattern, product) - before[name]))
                    for name, change in changes.items():
                        if name not in result["equations_by_reaction"]:
                            continue
                        row = result["equations_by_reaction"][name][reaction]
                        expected = sum(Fraction(c) * before[target] for target, c in row.items())
                        assert change == expected, (name, reaction, edges, opened, change, expected)
                        generator_checks += 1
    report = {"passed": True, "gated": gated, "feedback": feedback,
        "partial_equations_only": partial, "basis_size": len(forms), "states": states,
        "symbolic_generator_identities_all_n": symbolic,
        "independent_graph_pattern_counts": count_checks,
        "independent_graph_generator_identities": generator_checks,
        "independent_raw_propensity_checks": propensity_checks,
        "successor_graph_invariant_checks": products, "seconds": perf_counter() - start,
        "scope": "Symbolic oracle is all-n and independent of composition/counting. Full-graph checks are finite regressions. Invariant preservation and nonexplosion have separate mathematical arguments."}
    seed = forms[inputs.data["initial_observables"][0]]
    degree = int(sp.Poly(seed, n).degree())
    if degree >= 2 and gated and not feedback:
        report["independent_minimality_diagnostic"] = analytic_cyclic_rank(degree)
        # Exercise the actual nonzero high-order seed on larger complete hosts.
        # Counting ordered distinct isolated A-nodes is elementary label-count
        # combinatorics, not the footprint/partition implementation. We do not
        # pretend this is factorial-sized exhaustive match enumeration.
        seed_name = inputs.data["initial_observables"][0]
        coefficient, pattern = observables[seed_name][0]
        assert len(observables[seed_name]) == 1 and coefficient == 1
        assert not edge_rows(rule_graph(pattern)) and not engine_no_edges(pattern) and not engine_conditions(pattern)
        assert len(node_rows(rule_graph(pattern))) == degree and set(node_labels(rule_graph(pattern)).values()) == {"A"}
        large_counts = large_derivatives = 0
        for edges in (degree - 1, degree, degree + 2):
            for opened in (False, True):
                graph = host(edges, 0, opened)
                actual = perm(sum(label == "A" for label in node_labels(graph).values()), degree)
                assert actual > 0 and actual == seed.subs({n: edges, m: 0, x: int(opened)})
                large_counts += 1
                for declaration in inputs.data["rules"]:
                    reaction = declaration["id"]
                    change = Fraction()
                    for product, multiplicity in transitions(rules[reaction], graph).items():
                        target_count = perm(sum(label == "A" for label in node_labels(product).values()), degree)
                        change += multiplicity * Fraction(declaration["per_match"]) * (target_count - actual)
                    row = result["equations_by_reaction"][seed_name][reaction]
                    derived = sum((sp.Rational(c) * forms[target] for target, c in row.items()), sp.Integer(0))
                    assert change == derived.subs({n: edges, m: 0, x: int(opened)})
                    large_derivatives += 1
        report["nonzero_seed_large_host_checks"] = {
            "path_vertices": [degree, degree + 1, degree + 3],
            "seed_counts": large_counts, "seed_generator_identities": large_derivatives,
            "method": "Count ordered distinct isolated A-nodes as math.perm(actual_host_label_count, degree); apply every literal graph-rule match independently",
            "scope": "Seed only; RHS uses basis polynomials separately checked symbolically. Not exhaustive high-order seed-match enumeration and not full high-order basis enumeration."}
    report["seconds"] = perf_counter() - start
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("--ungated", action="store_true")
    parser.add_argument("--feedback", action="store_true")
    parser.add_argument("--partial", action="store_true")
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = verify(load(args.model), json.loads(args.result.read_text(encoding="utf-8")), max_n=args.max_n,
                    gated=not args.ungated, feedback=args.feedback, partial=args.partial)
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
