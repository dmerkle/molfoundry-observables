"""Multiplicity-preserving conditional-DPO compositions and jump observables."""
from collections import defaultdict
from dataclasses import dataclass, field
import molfoundry as mf
from mfobs.core.composition import CompositionDiagnostics, group_compositions, rule_labels
from mfobs.core.forbidden import incidence_constraints


def grouped_common_compositions(first, second, *, connected=False,
                                maximum=False, include_empty=True):
    diagnostics, counts, representatives = group_compositions(
        first, second, connected=connected, maximum=maximum, include_empty=include_empty)
    return diagnostics.accepted_raw, counts, representatives


def grouped_restricted_common_compositions(first, second, *, forbidden_specification,
                                           connected=False, maximum=False, include_empty=True):
    # Compiled constraints prune partial overlaps that cannot satisfy the language.
    constraints = incidence_constraints(forbidden_specification, *rule_labels(first, second))
    return group_compositions(first, second, constraints=constraints,
        connected=connected, maximum=maximum, include_empty=include_empty,
        maximum_intermediate_graphs=forbidden_specification.get('maximumCandidateOverlaps'))


def dpo_jump_closure(rule, name='jump'):
    return rule.dpoJumpClosure(name)


def sqpo_jump_closure(rule, name='jump'):
    return rule.sqpoJumpClosure(name)


@dataclass
class CommutatorResult:
    """One reaction's commutator with one observable.

    ``coefficient`` follows the Posselt/TCS table convention
    ``count(reaction * observable) - count(observable * reaction)``.
    """

    left_raw: int
    right_raw: int
    composite_coefficients: dict[str, int]
    composite_representatives: dict[str, object]
    observable_coefficients: dict[str, int]
    observable_representatives: dict[str, object]
    left_restricted: CompositionDiagnostics | None = None
    right_restricted: CompositionDiagnostics | None = None
    composition_orders: dict = field(default_factory=dict)
    raw_composites: list = field(default_factory=list)


def commutator(
    observable,
    reaction,
    *,
    connected: bool = False,
    include_empty: bool = True,
    forbidden_specification: dict | None = None,
    constraints=None,
    rule_rejection=None,
    record_compositions: bool = False,
) -> CommutatorResult:
    """Compute and canonically collect both sides of the commutator."""

    raw_composites, orders = [], {}
    if forbidden_specification is not None:
        constraints = [*(constraints or ()),
                       *incidence_constraints(forbidden_specification, *rule_labels(observable, reaction))]

    def compose(first, second, order, sign):
        def record(rule, reason):
            raw_composites.append({"order": order, "coefficient": sign,
                                   "gml": rule.getGMLString(), "rejection_reason": reason})
        diagnostics, counts, representatives = group_compositions(
            first, second, constraints=constraints, connected=connected,
            include_empty=include_empty,
            maximum_intermediate_graphs=(forbidden_specification.get('maximumCandidateOverlaps')
                                         if forbidden_specification is not None else None),
            rule_rejection=rule_rejection, on_rule=record if record_compositions else None)
        orders[order] = {"public_raw": diagnostics.accepted_raw + diagnostics.forbidden_input_rejected,
                         "retained": diagnostics.accepted_raw, "rejected": diagnostics.forbidden_reasons}
        return diagnostics, counts, representatives

    left, left_counts_all, left_reps_all = compose(observable, reaction, "before", -1)
    right, right_counts_all, right_reps_all = compose(reaction, observable, "after", 1)
    left_raw, right_raw = left.accepted_raw, right.accepted_raw
    left_restricted = left if forbidden_specification is not None else None
    right_restricted = right if forbidden_specification is not None else None
    left_counts = left_counts_all
    right_counts = right_counts_all
    left_reps = left_reps_all
    right_reps = right_reps_all
    all_composites = set(left_counts) | set(right_counts)
    composite_coefficients = {
        key: right_counts[key] - left_counts[key]
        for key in all_composites
        if right_counts[key] != left_counts[key]
    }
    composite_representatives = left_reps | right_reps

    observable_coefficients: defaultdict[str, int] = defaultdict(int)
    observable_representatives: dict[str, object] = {}
    for index, (key, coefficient) in enumerate(
        sorted(composite_coefficients.items())
    ):
        jump = dpo_jump_closure(
            composite_representatives[key], name=f"jump_{index}"
        )
        jump_key = jump.canonKey
        observable_coefficients[jump_key] += coefficient
        observable_representatives.setdefault(jump_key, jump)

    observable_coefficients = defaultdict(
        int,
        {
            key: coefficient
            for key, coefficient in observable_coefficients.items()
            if coefficient
        },
    )
    return CommutatorResult(
        left_raw=left_raw,
        right_raw=right_raw,
        composite_coefficients=composite_coefficients,
        composite_representatives=composite_representatives,
        observable_coefficients=dict(observable_coefficients),
        observable_representatives=observable_representatives,
        left_restricted=left_restricted,
        right_restricted=right_restricted,
        composition_orders=orders,
        raw_composites=raw_composites,
    )
