"""Finite-tree oracle and direct MolFoundry checks for Rémy's rules."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from functools import cache
from typing import Mapping

from mfobs.core.engine import mf
from mfobs.models.remy.inputs import Tree, LEAF, tree_gml, insertion_rules
from mfobs.models.remy.events import derive_jump_classes

@cache
def trees(internal_nodes: int) -> tuple[Tree, ...]:
    """Finite Catalan enumerator used only by validation and visual samples."""

    if internal_nodes < 0:
        raise ValueError("internal node count must be nonnegative")
    if internal_nodes == 0:
        return (LEAF,)
    return tuple(
        (left, right)
        for left_size in range(internal_nodes)
        for left in trees(left_size)
        for right in trees(internal_nodes - 1 - left_size)
    )


@cache
def observable_values(tree: Tree) -> tuple[int, int, int, int]:
    """Evaluate E,P1,P2,P3 from their recursive pattern definitions."""

    if not tree:
        return (1, 0, 0, 0)
    left, right = tree
    lv = observable_values(left)
    rv = observable_values(right)
    p1 = 1 + lv[1] + rv[1]
    p2 = lv[2] + rv[2] + int(bool(left))
    p3 = lv[3] + rv[3] + int(bool(left) and bool(left[0]))
    return (2 * p1 + 1, p1, p2, p3)


@dataclass(frozen=True)
class InsertionEvent:
    target: Tree
    path: str
    side: str
    edge_type: str


@cache
def insertion_events(tree: Tree) -> tuple[InsertionEvent, ...]:
    """All Rémy events, retaining equal targets as distinct matches."""

    result = [
        InsertionEvent((tree, LEAF), "", "L", "I"),
        InsertionEvent((LEAF, tree), "", "R", "I"),
    ]
    if tree:
        left, right = tree
        result.extend(
            InsertionEvent((event.target, right), "L" + event.path, event.side,
                           "L" if not event.path else event.edge_type)
            for event in insertion_events(left)
        )
        result.extend(
            InsertionEvent((left, event.target), "R" + event.path, event.side,
                           "R" if not event.path else event.edge_type)
            for event in insertion_events(right)
        )
    return tuple(result)


def _propensity(event: Mapping, values: tuple[int, int, int, int]) -> int:
    return sum(
        int(coefficient) * value
        for coefficient, value in zip(event["propensity"], values)
    )


def validate_finite(maximum_source_internal_nodes: int = 8) -> dict:
    """Independent bounded oracle; not used by symbolic derivation."""

    derivation = derive_jump_classes()
    classes = {tuple(event["delta"]): event for event in derivation["classes"]}
    source_count = 0
    rewrite_events = 0
    for size in range(maximum_source_internal_nodes + 1):
        for tree in trees(size):
            source_count += 1
            before = observable_values(tree)
            actual: Counter = Counter()
            for event in insertion_events(tree):
                after = observable_values(event.target)
                delta = tuple(a - b for a, b in zip(after, before))
                actual[delta] += 1
                rewrite_events += 1
            expected = {
                delta: _propensity(event, before)
                for delta, event in classes.items()
            }
            if actual != Counter(expected):
                raise AssertionError(
                    f"finite event partition mismatch at size {size}: {tree}"
                )
            if sum(actual.values()) != 2 * before[0]:
                raise AssertionError("event multiplicity does not equal 2E")
    return {
        "status": "PASS",
        "scope": "bounded validation only; excluded from all-tree derivation",
        "maximum_source_internal_nodes": maximum_source_internal_nodes,
        "maximum_target_internal_nodes": maximum_source_internal_nodes + 1,
        "source_trees": source_count,
        "rewrite_events_counting_multiplicity": rewrite_events,
        "all_equal_targets_retained": True,
    }


def validate_molfoundry_direct(maximum_source_internal_nodes: int = 3) -> dict:
    """Compare the six actual GML rules to the recursive Rémy transitions."""

    rules = insertion_rules()
    source_count = 0
    public_hyperedges = 0
    injective_matches = 0
    for size in range(maximum_source_internal_nodes + 1):
        for index, tree in enumerate(trees(size)):
            source_count += 1
            source = mf.graphGMLString(
                tree_gml(tree, name=f"source_{size}_{index}"), add=False
            )
            expected = Counter(
                mf.graphGMLString(tree_gml(event.target), add=False).canonKey
                for event in insertion_events(tree)
            )
            actual: Counter = Counter()
            dg = mf.DG(graphDatabase=[source])
            with dg.build() as builder:
                for rule in rules.values():
                    injective_matches += rule.numMatches([source], mf.ApplyMode.UseAny)
                    for edge in builder.apply(
                        [source], rule, mode=mf.ApplyMode.UseAny, verbosity=0
                    ):
                        targets = list(edge.targets)
                        if len(targets) != 1:
                            raise AssertionError("Rémy rule must have one graph product")
                        actual[targets[0].graph.canonKey] += edge.numMatches(rule)
                        public_hyperedges += 1
            if actual != expected:
                raise AssertionError(
                    f"MolFoundry/direct transition mismatch at size {size}: {tree}"
                )
    expected_matches = sum(
        len(insertion_events(tree))
        for size in range(maximum_source_internal_nodes + 1)
        for tree in trees(size)
    )
    if injective_matches != expected_matches:
        raise AssertionError("actual GML LHS match count differs from 2E")
    return {
        "status": "PASS",
        "maximum_source_internal_nodes": maximum_source_internal_nodes,
        "source_trees": source_count,
        "injective_LHS_matches": injective_matches,
        "public_DG_distinct_hyperedges": public_hyperedges,
        "public_DG_note": "DG multiplicities read with HyperEdge.numMatches",
        "match_multiplicities_preserved": True,
        "direction_and_order_encoded_by_incidence_roles_and_I_L_R_types": True,
    }


from mfobs.core.typed_incidence import typed_incidence_commutator
from mfobs.models.remy.inputs import TYPED_INCIDENCE_SPECIFICATION, observable_rules


def raw_first_commutators() -> dict:
    """Compute raw MolFoundry composite expressions before jump closure."""

    rules = insertion_rules()
    observables = observable_rules()
    result = {}
    for observable_name, components in observables.items():
        coefficients: Counter = Counter()
        representatives = {}
        rows = []
        for component in components:
            for reaction_name, reaction in rules.items():
                value = typed_incidence_commutator(
                    component,
                    reaction,
                    specification=TYPED_INCIDENCE_SPECIFICATION,
                    connected=False,
                    include_empty=True,
                )
                coefficients.update(value.coefficients)
                representatives.update(value.representatives)
                rows.append(
                    {
                        "observable_component": component.name,
                        "reaction": reaction_name,
                        "left": value.left.record(),
                        "right": value.right.record(),
                    }
                )
        coefficients = Counter(
            {key: value for key, value in coefficients.items() if value}
        )
        result[observable_name] = {
            "coefficients": dict(coefficients),
            "representatives": {
                key: representatives[key] for key in coefficients
            },
            "ordered_composition_rows": rows,
        }
    return result
