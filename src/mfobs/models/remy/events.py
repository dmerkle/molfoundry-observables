"""All-tree insertion-context counts and the resulting complete jump model."""
from __future__ import annotations
from functools import cache
import sympy as sp
from mfobs.core.analysis_contract import AnalysisContract, NormalizedRule, ObservableDeclaration, StateLanguageContract
from mfobs.core.composition_jump_model import (
    CompositionJumpCompilation,
    GraphReactionClass,
    compile_composition_jump_model,
)
from mfobs.core.jump_closure_worklist import RuleExpression
from mfobs.core.polynomial_jump_closure import JumpEvent, JumpModelProvenance, PolynomialJumpModel, PolynomialRing
from mfobs.core.typed_incidence import pair_commutator
from mfobs.models.remy.inputs import EDGE_TYPES, TYPED_INCIDENCE_SPECIFICATION, insertion_rule_gml, insertion_rules, observable_names, observable_rules
from mfobs.models.remy.counting import REMY_JUMP_PROJECTION, _left_spine_component, valid_tree_dpo_jump_closure


@cache
def composition_prefix_compilation(order: int = 3) -> CompositionJumpCompilation:
    """Derive the universally constant ``(E,P1)`` jump prefix by DPO and valid-tree count reduction.

    The original six rules are already constant-jump classes for these two
    coordinates.  Higher ``P_k`` jumps depend on tree context and therefore
    remain the responsibility of the separately certified all-tree partition.
    """

    rules = insertion_rules()
    all_observables = observable_rules()
    coordinates = {
        name: RuleExpression.from_rules(all_observables[name])
        for name in ("E", "P1")
    }
    reactions = tuple(
        GraphReactionClass(
            identifier=name,
            channel="rho",
            rule=rule,
            source_rule=name,
            predicate=f"selected {name[-1]} edge; old subtree on {name[1]}",
            metadata={"side": name[1], "stem_type": name[-1]},
        )
        for name, rule in rules.items()
    )

    def activity(ring, reaction, activity_rule):
        depth, stem = _left_spine_component(activity_rule)
        if depth != 0 or stem != reaction.metadata["stem_type"]:
            raise ValueError("unexpected reduced DPO Rémy activity identity")
        return (
            ring.constant(1)
            if stem == "I"
            else ring.coordinate("P1")
        )

    return compile_composition_jump_model(
        name="remy-tree-E-P1-prefix",
        ring=PolynomialRing(("E", "P1")),
        observables=coordinates,
        reactions=reactions,
        pair_commutator=pair_commutator(TYPED_INCIDENCE_SPECIFICATION, include_empty=False),
        activity_reducer=activity,
        event_partition_certificate={
            "status": "VERIFIED_ALL_STATES",
            "method": "six typed edge/side rule classes",
            "argument": (
                "Every valid tree edge has exactly one I/L/R type and each "
                "selected edge has exactly the two GL/GR insertion choices."
            ),
            "enumeration_used": False,
        },
        input_definition_sha256=analysis_contract(order).digest,
        jump_closure=valid_tree_dpo_jump_closure,
        jump_projection_name=REMY_JUMP_PROJECTION,
        composition_backend_name=(
            "MolFoundry DPO typed-incidence composition with valid-tree count reduction"
        ),
    )


def _local_jump(side: str, below_depth: int, left_ancestors: int, order: int) -> tuple[int, ...]:
    inserted_depth = 1 + below_depth if side == "L" else 1
    result = []
    for threshold in range(1, order + 1):
        change = int(inserted_depth >= threshold)
        change += sum(
            int(distance + inserted_depth >= threshold)
            - int(distance + below_depth >= threshold)
            for distance in range(1, left_ancestors + 1)
        )
        result.append(change)
    return tuple(result)


def derive_jump_classes(order: int = 3) -> dict:
    """Derive event propensities symbolically for all valid binary trees.

    Let ``d`` be the left-spine depth below the selected edge and ``a`` the
    consecutive run of left-child edges above it.  The all-tree bijection
    ``#{edges | d>=r and a>=s} = P_(r+s)`` supplies every category count.
    Möbius inversion then partitions events without enumerating trees.
    """

    if order < 1:
        raise ValueError("at least P1 is required")
    cap = order + 2
    symbols = sp.symbols(f"P0:{2 * cap + 2}")

    def survival(below: int, above: int):
        return symbols[below + above]

    def exact_bin(below: int, above: int):
        value = survival(below, above)
        if below < cap:
            value -= survival(below + 1, above)
        if above < cap:
            value -= survival(below, above + 1)
        if below < cap and above < cap:
            value += survival(below + 1, above + 1)
        return sp.expand(value)

    grouped: dict[tuple[int, ...], sp.Expr] = {}
    categories = []
    for below in range(cap + 1):
        for above in range(cap + 1):
            mass = exact_bin(below, above)
            for side in ("L", "R"):
                delta = _local_jump(side, below, above, order)
                grouped[delta] = sp.expand(grouped.get(delta, 0) + mass)
                categories.append(
                    {
                        "below_bin": f">={cap}" if below == cap else str(below),
                        "above_bin": f">={cap}" if above == cap else str(above),
                        "side": side,
                        "delta_P": list(delta),
                        "symbolic_mass": str(mass),
                    }
                )

    classes = []
    unresolved = []
    for index, (delta_p, expression) in enumerate(grouped.items()):
        expression = sp.expand(expression)
        polynomial = sp.Poly(expression, *symbols)
        coefficients = [polynomial.coeff_monomial(symbol) for symbol in symbols]
        residual = sp.expand(
            expression
            - sum(coefficient * symbol for coefficient, symbol in zip(coefficients, symbols))
        )
        if residual != 0 or any(coefficients[order + 1:]):
            unresolved.append({
                "delta": [2, *delta_p],
                "propensity_symbolic": str(expression).replace("P0", "E"),
                "required_observables": [
                    str(symbols[i]) for i in range(order + 1, len(symbols))
                    if coefficients[i]
                ],
                "nonlinear_residual": str(residual),
            })
            continue
        values = [int(coefficients[i]) for i in range(order + 1)]
        classes.append(
            {
                "id": f"J{index}",
                "delta": [2, *delta_p],
                "propensity": values,
                "propensity_symbolic": str(expression).replace("P0", "E"),
            }
        )

    total = sp.expand(sum(grouped.values()))
    if total != 2 * symbols[0]:
        raise AssertionError("symbolic event partition does not sum to 2E")
    record = {
        "status": "VERIFIED_ALL_VALID_TREES",
        "method": "finite-context partition plus exact symbolic Mobius inversion",
        "enumeration_used": False,
        "observable_family": {
            "P_k": (
                "number of internal vertices whose rooted subtree has a "
                "left spine containing at least k internal vertices"
            ),
            "P0": "E, the number of selectable directed typed edges",
        },
        "universal_bijection": "N(r,s) = #{selected edges: d>=r and a>=s} = P_(r+s)",
        "bijection_proof": {
            "forward": (
                "for s>0 map a selected edge to its ancestor s left edges above; "
                "for s=0 map it to its source subtree root"
            ),
            "inverse": (
                "from a P_(r+s) occurrence descend s ordered L edges to the "
                "unique selected edge; the two maps are inverse"
            ),
            "boundary_case": "N(0,0)=E",
        },
        "proof_boundary": (
            "the tree-language bijection is model-specific semantic input; "
            f"the observable basis itself is discovered from P{order} by the generic "
            "first-commutator/jump-projection worklist; Mobius inversion, cancellation, span "
            "closure, and EMGF export are machine-checked"
        ),
        "local_jump_formula": (
            "Delta P_k = [d' >= k] + sum_{j=1..a} "
            "([j+d' >= k]-[j+d >= k]), with d'=d+1 for GL and d'=1 for GR"
        ),
        "partition_cap": cap,
        "higher_observable_coefficients_cancel": True,
        "total_event_propensity": "2*E",
        "classes": classes,
        "category_derivation": categories,
    }
    if unresolved:
        record.update(
            status="UNKNOWN_UNSUPPORTED_EVENT_REDUCTION",
            reason="Event activities require counts outside the selected left-spine coordinates.",
            higher_observable_coefficients_cancel=False,
            unresolved=unresolved,
        )
        # Partial classes do not form an event model on these coordinates.
        record.pop("classes")
    return record


@cache
def analysis_contract(order: int = 3) -> AnalysisContract:
    """Normalized all-tree input shared by linear discovery and pJC/EMGF."""

    rules = tuple(
        NormalizedRule(
            identifier=f"G{side}_{stem}",
            rate="rho",
            semantics="conditional-DPO",
            gml=insertion_rule_gml(side, stem),
        )
        for side in ("L", "R")
        for stem in EDGE_TYPES
    )
    names = observable_names(order)
    selected = f"P{order}"
    observables = tuple(
        ObservableDeclaration(
            name,
            (
                "number of selectable directed typed edges"
                if name == "E"
                else f"number of ordered full-binary left-spine patterns of depth {name[1:]}"
            ),
            "selected" if name == selected else "derived",
            (
                "literal three-component identity-rule expression"
                if name == selected
                else "first-commutator residual, structurally named after discovery"
            ),
        )
        for name in names
    )
    language = StateLanguageContract(
        identifier="directed-planar-full-binary-trees-v1",
        carrier=(
            "all finite rooted planar full binary trees with one exterior I stem "
            "and directed typed incidence edges I/L/R"
        ),
        infinite=True,
        invariants=(
            "one exterior I stem and one root",
            "every internal vertex has exactly one ordered L and one ordered R child",
            "leaves have no children",
            "all child edges point toward the root and the directed graph is acyclic",
        ),
        overlap_constraints=(
            "each typed edge-node has exactly one src and one tgt incidence",
            "at most one outgoing represented edge per object",
            "at most one incoming edge of each type I/L/R",
            "simple directed and underlying graphs; no successor beyond an I edge",
            "every intermediate overlap and composite side satisfies the tree theory",
        ),
        preservation_argument=(
            "Each of the six direct rules replaces one selected edge by one new "
            "binary vertex, one new leaf and three typed edges. The structural "
            "overlap filter preserves functionality, order, the unique stem and acyclicity."
        ),
        semantic_quotient=(
            "event classes induced by below-left-spine depth and consecutive "
            "left ancestors, reduced by N(r,s)=P_(r+s)"
        ),
    )
    return AnalysisContract(
        model_id="remy-trees",
        title="Rémy trees from six typed insertion rules",
        rules=rules,
        state_language=language,
        observables=observables,
        linear_seeds=(selected,),
        pjc_targets=names,
        normalization={
            "rule_matches": "one event per injective selected-edge/side match",
            "typed_edges": "undirected incidence encoding preserves source and target ports",
            "rewrite_semantics": "conditional DPO on valid incidence trees; Behr's SqPO insertions have the same matches and targets here",
            "jump_projection": REMY_JUMP_PROJECTION,
        },
    )


@cache
def polynomial_jump_model(order: int = 3) -> PolynomialJumpModel:
    """Convert the symbolically derived all-tree event partition to generic pJC."""

    contract = analysis_contract(order)
    names = observable_names(order)
    ring = PolynomialRing(names)
    events = []
    derivation = derive_jump_classes(order)
    if derivation["status"] != "VERIFIED_ALL_VALID_TREES":
        from mfobs.core.pattern_support import UnsupportedReduction
        raise UnsupportedReduction(derivation["reason"])
    for event in derivation["classes"]:
        propensity = ring.add(*(
            ring.scale(ring.coordinate(name), coefficient)
            for name, coefficient in zip(names, event["propensity"])
            if coefficient
        ))
        events.append(JumpEvent(
            identifier=event["id"],
            channel="rho",
            delta=tuple(event["delta"]),
            propensity=propensity,
            rule="sum of G{L/R}_{I/L/R} matches in this event class",
            predicate="symbolic (below-depth, left-ancestor) event class",
        ))
    certificate = {
        "status": "VERIFIED_ALL_STATES",
        "method": "finite-context partition plus exact symbolic Mobius inversion",
        "argument": (
            "The bijection N(r,s)=P_(r+s) partitions every selected edge and side "
            "of every valid finite tree. Exact class propensities sum to 2*E; "
            "no enumerated tree is used to construct the partition."
        ),
        "enumeration_used": False,
    }
    return PolynomialJumpModel(
        contract.model_id,
        ring,
        tuple(events),
        certificate,
        contract.digest,
        JumpModelProvenance(
            backend="symbolic all-tree event quotient after rule-composition discovery",
            event_partition_origin=(
                "finite local (d,a) partition plus exact all-tree Mobius inversion"
            ),
            jump_origin="local tree insertion evaluated on generated P_k coordinates",
            activity_origin="all-tree bijection N(r,s)=P_(r+s)",
            jump_projection=REMY_JUMP_PROJECTION,
            composition_scope=(
                f"E/P1 jumps and P{order}-only observable discovery; "
                "higher context-dependent event jumps use the all-tree theorem"
            ),
            all_jumps_composition_derived=False,
            supplied_delta_table=False,
        ),
    )
