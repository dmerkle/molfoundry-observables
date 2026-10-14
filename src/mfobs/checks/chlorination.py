"""Direct MolFoundry applications and an independent four-site reference check."""
from __future__ import annotations
from collections import Counter, defaultdict

import molfoundry as mf
import sympy as sp


def _enabled_positions(mask: int):
    """Allowed neutral-ring positions for the literal NAC rule."""

    result = []
    for index in range(4):
        if mask & (1 << index):
            continue
        if index > 0 and mask & (1 << (index - 1)):
            continue
        if index < 3 and mask & (1 << (index + 1)):
            continue
        result.append(index)
    return result


def _molecule(mask: int, n: int = 0):
    """Concrete member of the unbounded-tail family, for post-checks only."""
    from mfobs.models.chlorination.inputs import molecule
    return molecule(n=n, mask=mask)


def _hcl():
    return mf.Graph.fromGMLString(
        'graph [ node [ id 0 label "H" ] node [ id 1 label "Cl" ] '
        'edge [ source 0 target 1 label "-" ] ]',
        add=False,
    )


def validate_all_masks(backend, result):
    """Finite differential checks; the all-tail claim rests on locality."""
    from mfobs.models.chlorination.inputs import methane as methane_graph, h2 as hydrogen_graph

    from types import SimpleNamespace
    if result["pjc"]["status"] != "CLOSED_FINITE_EVENT_MODEL":
        return {"status": "SKIPPED", "reason": result["pjc"].get("reason") or "incomplete synthesis"}
    names = result["pjc"]["coordinates"]
    symbols = {n: sp.Symbol(n) for n in names}
    engine_images = {n: sp.sympify(row["chlorinate"]["text"], locals=symbols)
                     for n, row in result["synthesis"]["support"]["generator_images"].items()}
    events = [SimpleNamespace(activity=sp.sympify(e["activity"]["text"], locals=symbols),
                              delta=tuple(e["delta"].get(n, 0) for n in names))
              for e in result["events"]]

    failures = []
    scalar_checks = reachability_checks = growth_checks = 0
    names = tuple(engine_images)
    patterns = {n: backend.patterns[n] for n in names}
    # n counts carbons beyond the methyl stub. These are regression samples,
    # not a bound on the language or evidence by extrapolation.
    contexts = ((0, 0, 0), (1, 2, 3), (7, 3, 1))
    methane = methane_graph()
    hydrogen = hydrogen_graph()

    def concrete_counts(graph):
        # Positive observables are identity rules. MolFoundry counts their
        # actual embeddings, including automorphisms, through the public API.
        counts = {}
        for name, pattern in patterns.items():
            dg = mf.DG(graphDatabase=[graph])
            with dg.build() as builder:
                edges = builder.apply([graph], pattern.rule, mode=mf.ApplyMode.UseAny, verbosity=0)
            counts[name] = sum(edge.numMatches(pattern.rule) for edge in edges)
        return counts

    def apply_heads(components, rule, molecule, reagent, inventory, before, factor=1):
        dg = mf.DG(graphDatabase=components)
        with dg.build() as builder:
            derivations = builder.apply(components, rule, mode=mf.ApplyMode.UseAny, verbosity=0)
        heads, jumps = set(), defaultdict(int)
        for edge in derivations:
            # Both rules consume one head and one distinct reagent species.
            assert Counter(v.graph.canonKey for v in edge.sources) == Counter(
                [molecule.canonKey, reagent.canonKey])
            after = dict.fromkeys(names, 0)
            for vertex in edge.targets:
                graph = vertex.graph
                if any(atom.stringLabel == "C-" for atom in graph.vertices):
                    heads.add(graph.canonKey)
                for name, value in concrete_counts(graph).items():
                    after[name] += value
            delta = tuple(after[name] - before[name] for name in names)
            jumps[delta] += edge.numMatches(rule) * inventory * factor
        return heads, dict(jumps)

    for n, methane_count, h2_count in contexts:
        values_by_mask = {
            mask: _coordinate_values(mask, patterns, n=n, methane=methane_count,
                                     h2=h2_count, hcl=backend.chlorine_total-mask.bit_count())
            for mask in range(16)}
        for mask in range(16):
            values = values_by_mask[mask]
            hcl = backend.chlorine_total - values["Z"]
            enabled = (_enabled_positions(mask) if backend.guarded else
                       [i for i in range(4) if not mask & (1 << i)])
            direct = {name: 0 for name in names}
            direct_events = defaultdict(int)
            for index in enabled:
                after = values_by_mask[mask | (1 << index)]
                delta = tuple(after[name] - values[name] for name in names)
                direct_events[delta] += hcl
                for name in direct:
                    direct[name] += hcl * (after[name] - values[name])

            for name, expected in direct.items():
                actual = _eval(engine_images[name], values)
                scalar_checks += 1
                if actual != expected:
                    failures.append({"mask": mask, "n": n, "kind": "generator",
                                     "coordinate": name, "expected": expected, "actual": actual})
            model_events = defaultdict(int)
            for event in events:
                model_events[event.delta] += _eval(event.activity, values)
            model_events = {delta: rate for delta, rate in model_events.items() if rate}
            direct_events = {delta: rate for delta, rate in direct_events.items() if rate}
            if model_events != direct_events:
                failures.append({"mask": mask, "n": n, "kind": "event-partition",
                                 "expected": {str(k): v for k, v in direct_events.items()},
                                 "actual": {str(k): v for k, v in model_events.items()}})

            molecule = _molecule(mask, n=n)
            # Repeated molecules have equal graph keys. Reachability needs one
            # representative of each available component type; rates above use
            # their actual inventories; numMatches supplies the per-species multiplicities.
            components = [molecule] + ([_hcl()] if hcl else [])
            components += ([methane] if methane_count else []) + ([hydrogen] if h2_count else [])
            before = concrete_counts(molecule)
            if before != values:
                failures.append({"mask": mask, "n": n, "kind": "molfoundry-counts",
                                 "expected": values, "actual": before})
            actual_heads, actual_events = apply_heads(
                components, backend.rule, molecule, _hcl(), hcl, before)
            if actual_events != model_events:
                failures.append({"mask": mask, "n": n, "kind": "molfoundry-joint-rates",
                                 "expected": {str(k): str(v) for k, v in model_events.items()},
                                 "actual": {str(k): str(v) for k, v in actual_events.items()}})
            expected_heads = {_molecule(mask | (1 << i), n=n).canonKey for i in enabled} if hcl else set()
            reachability_checks += 1
            if actual_heads != expected_heads:
                failures.append({"mask": mask, "n": n, "kind": "molfoundry-reachability",
                                 "expected": sorted(expected_heads), "actual": sorted(actual_heads)})

            grown_heads, growth_events = apply_heads(
                components, backend.growth_rule, molecule, methane, methane_count,
                before, backend.growth_rate_factor)
            expected_events = {(0,) * len(names): methane_count} if methane_count else {}
            if growth_events != expected_events:
                failures.append({"mask": mask, "n": n, "kind": "molfoundry-growth-rates",
                                 "actual": {str(k): str(v) for k, v in growth_events.items()}})
            expected_growth = {_molecule(mask, n=n+1).canonKey} if methane_count else set()
            growth_checks += 1
            if grown_heads != expected_growth:
                failures.append({"mask": mask, "n": n, "kind": "growth-reachability",
                                 "expected": sorted(expected_growth), "actual": sorted(grown_heads)})
            if methane_count:
                after_growth = _coordinate_values(mask, patterns, n=n+1,
                    methane=methane_count-1, h2=h2_count+1, hcl=hcl)
                if values != after_growth:
                    failures.append({"mask": mask, "n": n, "kind": "growth-changed-observable"})

    return {
        "status": "PASS" if not failures else "FAIL",
        "role": "post-synthesis validation only; not used to construct support, equations, or events",
        "masks_checked": 16, "states_checked": 16 * len(contexts),
        "contexts": [{"n": n, "methane": m, "h2": p} for n, m, p in contexts],
        "generator_scalar_checks": scalar_checks,
        "molfoundry_joint_rate_checks": reachability_checks,
        "molfoundry_growth_rate_checks": growth_checks,
        "molfoundry_reachability_checks": reachability_checks,
        "growth_reachability_checks": growth_checks,
        "failures": failures,
    }


from mfobs.numerics.context_nac import _coordinate_values


def _eval(expression, values):
    value = expression.subs({sp.Symbol(n): v for n, v in values.items()})
    if not value.is_Integer:
        raise AssertionError(f"Nonintegral concrete count: {value}")
    return int(value)
