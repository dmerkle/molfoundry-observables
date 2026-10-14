"""Shared-contract analysis for all two-regular eight-slot topologies.

The existing composition worklist discovers 48 coordinates of semantic rank
47 from the sole occupied-slot query.  This adapter removes the one verified
post-synthesis dependency, evaluates the resulting independent coordinate
lift on all 83 symmetry classes, and delegates literal rule application plus
pJC/EMGF derivation to the same generic finite-language backends as the fixed ring.
"""
from __future__ import annotations

from mfobs.core.graph import node_rows, make_graph, edge_rows, has_edge, node_labels
from mfobs.core.engine import mf
from fractions import Fraction
from itertools import combinations_with_replacement, product
from pathlib import Path
from typing import Mapping

from mfobs.core.analysis_contract import (
    AnalysisContract,
    NormalizedRule,
    ObservableDeclaration,
    StateLanguageContract,
    assert_same_input,
)
from mfobs.core.conditions import matches
from mfobs.core.finite_language_pjc import compile_finite_language_jump_model
from mfobs.core.finite_rule_projection import (
    FiniteRuleChannel,
    project_finite_rule_generators,
)

from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
from mfobs.core.gml import read_rule_record
from mfobs.inputs import load
from mfobs.models import ring_language as ring


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MODEL_PATH = PROJECT_ROOT / "examples" / "eight-slots" / "model.json"
TOPOLOGIES = ((8,), (4, 4), (3, 5))


def analysis_contract() -> AnalysisContract:
    inputs = load(MODEL_PATH)
    invariants = inputs.read("state_language")
    observable = inputs.data["observables"][0]
    return AnalysisContract(
        model_id="eight-slots-all-topologies",
        title="Eight locally owned slots on every two-regular topology",
        rules=tuple(
            NormalizedRule(
                record["id"], record["rate"], "conditional-DPO",
                inputs.file(record["file"]).read_text(encoding="utf-8"),
                (("NoEdge forbids a pre-existing bound edge",)
                 if record["id"] != "unbind" else ()),
            )
            for record in inputs.data["rules"]
        ),
        state_language=StateLanguageContract(
            identifier=invariants["kind"],
            carrier=(
                "all finite simple two-regular graphs of exactly eight owned C/P/L "
                "slots; connectedness is not required"
            ),
            infinite=False,
            invariants=(
                "exactly eight C centres, each owning one P and one L",
                "every C has exactly two distinct next neighbours",
                "each owned P/L pair carries at most one bound edge",
                "simple graph; no components outside the owned-slot language",
            ),
            overlap_constraints=(
                "all common overlaps, including empty and disconnected overlaps",
                "injective typed matches and transported NoEdge conditions",
                "reject ownership violations, branching and incomplete slots",
                "preserve every distinct neighbour match as stochastic multiplicity",
            ),
            preservation_argument=(
                "All three rules preserve every C/P/L and next edge and only toggle "
                "one owned bound edge; therefore component topology is invariant."
            ),
            semantic_quotient=(
                "83 complete occupancy classes: 30 for 8, 21 for 4+4 and 32 for 3+5"
            ),
            finite_enumeration_role=(
                "complete proof domain for this finite state language; discovery itself "
                "uses symbolic component-polynomial reduction"
            ),
        ),
        observables=(ObservableDeclaration(
            "P0", "number of occupied slots", "selected",
            "sole checked-in identity-rule query",
            inputs.file(observable["file"]).read_text(encoding="utf-8"),
        ),),
        linear_seeds=("P0",),
        pjc_targets=("P0",),
        normalization={
            **{record["id"]: f'{record["per_match"]} per injective match'
               for record in inputs.data["rules"]},
            "coordinate_dependency": (
                "one exact dependency is removed only after the 48-coordinate "
                "composition worklist closes"
            ),
        },
    )


def _cycle_classes(length: int) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted({
        ring.canonical_orbit(tuple(bits))
        for bits in product((0, 1), repeat=length)
    }))


def quotient_states() -> tuple[tuple[tuple[int, ...], ...], ...]:
    cycles = {length: _cycle_classes(length) for length in (3, 4, 5, 8)}
    states = [(value,) for value in cycles[8]]
    states.extend(tuple(pair) for pair in combinations_with_replacement(cycles[4], 2))
    states.extend((left, right) for left in cycles[3] for right in cycles[5])
    result = tuple(tuple(tuple(cycle) for cycle in state) for state in states)
    if len(result) != 83 or len(set(result)) != 83:
        raise AssertionError("eight-slot quotient must have 83 unique states")
    return result


def _representative(state: tuple[tuple[int, ...], ...]) -> mf.Graph:
    nodes = {}
    edges = []
    offset = 0
    for cycle in state:
        for local, occupied in enumerate(cycle):
            site = offset + local
            center, port, ligand = 3 * site, 3 * site + 1, 3 * site + 2
            nodes.update({center: "C", port: "P", ligand: "L"})
            edges.extend(((center, port, "cp"), (center, ligand, "cl")))
            if occupied:
                edges.append((port, ligand, "bound"))
        for local in range(len(cycle)):
            left = offset + local
            right = offset + (local + 1) % len(cycle)
            edges.append((3 * left, 3 * right, "next"))
        offset += len(cycle)
    if offset != 8:
        raise ValueError("quotient representative must contain eight slots")
    return make_graph(nodes, edges)


def _project(graph: mf.Graph) -> tuple[tuple[int, ...], ...]:
    labels = node_labels(graph)
    centers = {node for node, label in node_rows(graph) if label == "C"}
    nexts = {center: set() for center in centers}
    occupied = {}
    for center in centers:
        adjacency = [
            (right if left == center else left, label)
            for left, right, label in edge_rows(graph) if center in (left, right)
        ]
        ports = [node for node, label in adjacency
                 if label == "cp" and labels[node] == "P"]
        ligands = [node for node, label in adjacency
                   if label == "cl" and labels[node] == "L"]
        if len(ports) != 1 or len(ligands) != 1:
            raise ValueError("generated graph violates slot ownership")
        occupied[center] = int(has_edge(graph, ports[0], ligands[0], "bound"))
        nexts[center] = {
            node for node, label in adjacency
            if label == "next" and labels[node] == "C"
        }
        if len(nexts[center]) != 2:
            raise ValueError("generated graph violates two-regularity")

    unseen = set(centers)
    cycles = []
    while unseen:
        start = min(unseen)
        previous = None
        current = start
        ordered = []
        while current not in ordered:
            ordered.append(current)
            following = min(node for node in nexts[current] if node != previous)
            previous, current = current, following
        if current != start:
            raise ValueError("next component is not a simple cycle")
        unseen.difference_update(ordered)
        cycles.append(ring.canonical_orbit(tuple(occupied[node] for node in ordered)))
    return tuple(sorted(cycles, key=lambda cycle: (len(cycle), cycle)))


def _observable_vectors(data: Mapping, states) -> dict[str, tuple[Fraction, ...]]:
    vectors = {}
    hosts = tuple(_representative(state) for state in states)
    for name in data.get("basis", data["observable_expressions"]):
        terms = data["observable_expressions"][name]
        values = [Fraction()] * len(hosts)
        for term in terms:
            query = read_rule_record(term["gml"])
            values = [value + Fraction(term["coefficient"]) * len(
                matches(query, host))
                for value, host in zip(values, hosts)]
        vectors[name] = tuple(values)
    return vectors


def _coordinate_basis(vectors):
    import sympy as sp
    names = tuple(vectors)
    reduced, pivots = sp.Matrix([vectors[name] for name in names]).T.rref()
    dependencies = []
    for free in set(range(len(names))) - set(pivots):
        relation = {names[free]: "1"}
        relation.update({names[pivot]: str(-reduced[row, free])
                         for row, pivot in enumerate(pivots) if reduced[row, free]})
        dependencies.append(relation)
    return tuple(names[index] for index in pivots), dependencies


def finite_rule_projection(states=None):
    inputs = load(MODEL_PATH)
    states = quotient_states() if states is None else tuple(states)
    return project_finite_rule_generators(
        states=states,
        channels=tuple(
            FiniteRuleChannel(
                record["id"], read_rule_record(inputs.file(record["file"]).read_text(encoding="utf-8")),
                Fraction(record["per_match"]),
            )
            for record in inputs.data["rules"]
        ),
        representative=lambda state, _channel: _representative(state),
        project_target=lambda _source, _channel, target: _project(target),
    )


def _validate_first_moments(names, all_vectors, generators, result) -> dict:
    rows = {(name, channel): {target: Fraction(value) for target, value in coefficients.items()}
            for name, channels in result["equations_by_reaction"].items()
            for channel, coefficients in channels.items()}
    state_count = len(next(iter(all_vectors.values())))
    checks = 0
    for observable in names:
        for channel, generator in generators.items():
            expected = tuple(
                sum((generator[source][target] * all_vectors[observable][target]
                     for target in range(state_count)), Fraction())
                for source in range(state_count)
            )
            coefficients = rows.get((observable, channel), {})
            actual = tuple(
                sum((coefficient * all_vectors[target][source]
                     for target, coefficient in coefficients.items()), Fraction())
                for source in range(state_count)
            )
            if actual != expected:
                raise AssertionError(f"eight-slot first-moment mismatch: {observable}/{channel}")
            checks += state_count
    return {"status": "PASS", "exact_scalar_checks": checks}


def analyze(*, composition_result: Mapping) -> dict:
    """Derive finite events from an explicit closed composition result."""
    contract = analysis_contract()
    data = dict(composition_result)
    inputs = load(MODEL_PATH)
    if data.get("input_definition_sha256") != inputs.definition_digest():
        raise ValueError("eight-site composition result does not match its manifest")
    if not data.get("closed") or data.get("pending"):
        raise ValueError("eight-site composition worklist is not exactly closed")
    states = quotient_states()
    all_vectors = _observable_vectors(data, states)
    names, dependencies = _coordinate_basis(all_vectors)
    if len(data["basis"]) != 48 or len(names) != 47:
        raise ValueError("expected the 48-coordinate rank-47 eight-site result")
    projection = finite_rule_projection(states)
    compilation = compile_finite_language_jump_model(
        name=contract.model_id,
        states=states,
        observable_names=names,
        observable_vectors=tuple(all_vectors[name] for name in names),
        generators=projection.generators,
        input_definition_sha256=contract.digest,
        quotient_certificate={
            "method": "complete 83-class all-topology quotient",
            "topology_sector_sizes": {"8": 30, "4+4": 21, "3+5": 32},
            "finite_enumeration_role": contract.state_language.finite_enumeration_role,
        },
        source_rules={channel: channel for channel in projection.generators},
    )
    pjc = analyze_polynomial_jump_closure(compilation.model)
    pjc.update({
        "scope": "exact finite-language theorem on the complete eight-slot quotient",
        "selected_targets": list(contract.pjc_targets),
        "coordinate_lift": {
            "selected": list(contract.pjc_targets),
            "generated_support": [name for name in names if name != "P0"],
            "dropped_redundant_coordinate": next(
                name for name in data["basis"] if name not in names
            ),
            "origin": "closed composition worklist followed by its verified rank test",
        },
        "finite_language_compilation": compilation.evidence,
    })
    linear = {
        "schema_version": 1,
        "status": "CLOSED",
        "analysis": "linear-observable-discovery",
        "input_definition_sha256": contract.digest,
        "seeds": ["P0"],
        "basis": list(data["basis"]),
        "dimension": 48,
        "semantic_rank": 47,
        "scope": "exact modulo the supplied eight-slot ownership axioms",
        "ordered_composition_calls": data["ordered_composition_calls"],
        "raw_compositions": data["raw_compositions"],
    }
    assert_same_input(contract, (linear, pjc))
    return {
        "schema_version": 1,
        "contract": contract.record(),
        "linear": linear,
        "pjc": pjc,
        "finite_rule_projection": projection.evidence,
        "cross_branch_validation": _validate_first_moments(
            names, all_vectors, projection.generators, data
        ),
        "post_synthesis_dependency": dependencies[0],
    }
