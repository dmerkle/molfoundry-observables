"""Assemble linear and event results for the Rémy model."""
from __future__ import annotations
from typing import Mapping
from mfobs.core.analysis_contract import assert_same_input
from mfobs.core.operator_span import (
    discover_diagonal_commutator_span,
    emgf_operator_record,
)
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
from mfobs.models.remy.inputs import OBSERVABLES, OMEGA, VARIABLES
from mfobs.models.remy.counting import REMY_JUMP_PROJECTION, observable_synthesis
from mfobs.models.remy.events import analysis_contract, composition_prefix_compilation, derive_jump_classes, polynomial_jump_model


def common_analysis(*, linear_limit: int = 16, order: int = 3, record_patterns: bool = False, state_space=None) -> dict:
    """Run the common two-branch architecture on the Rémy model."""
    if state_space is not None:
        state_space.require(OMEGA)
    contract = analysis_contract(order)
    synthesis = observable_synthesis(order, maximum_dimension=linear_limit)
    synthesis_record = synthesis.record()
    linear = {
        "schema_version": 1,
        "status": synthesis.status,
        "analysis": "linear-observable-discovery",
        "algorithm": synthesis_record["algorithm"],
        "jump_projection": synthesis_record["jump_projection"],
        "input_definition_sha256": contract.digest,
        "seeds": list(contract.linear_seeds),
        "basis": list(synthesis.basis),
        "origins": synthesis.origins,
        "steps": synthesis.steps,
        "pending": synthesis.pending,
        "obligations": synthesis.obligations,
        "dimension": len(synthesis.basis),
        "raw_rule_expressions": {
            name: {**expression.record(), **({"gml_terms": [
                {"coefficient": str(expression.coefficients[key]),
                 "gml": expression.representatives[key].getGMLString()}
                for key in sorted(expression.coefficients)
            ]} if record_patterns else {})}
            for name, expression in synthesis.basis.items()
        },
        "scope": "exact on the declared all-tree typed-incidence language",
    }
    prefix = composition_prefix_compilation(order)
    derivation = derive_jump_classes(order)
    if derivation["status"] == "VERIFIED_ALL_VALID_TREES":
        pjc = analyze_polynomial_jump_closure(polynomial_jump_model(order))
    else:
        pjc = {
            "schema_version": 1,
            "analysis": "polynomial-jump-closure",
            "input_definition_sha256": contract.digest,
            "status": derivation["status"],
            "reason": derivation["reason"],
            "unresolved": derivation["unresolved"],
        }
    assert_same_input(contract, (linear, pjc))
    return {
        "contract": contract.record(),
        "composition": {
            **prefix.evidence,
            "scope": (
                "universally constant E/P1 jump prefix; higher P_k require the "
                "all-tree context partition recorded by the pJC certificate"
            ),
            "full_event_partition_composition_derived": False,
        },
        "model_specific_event_reduction": {
            key: derivation[key] for key in (
                "status", "method", "universal_bijection", "bijection_proof",
                "local_jump_formula", "proof_boundary", "reason", "unresolved",
            ) if key in derivation
        },
        "linear": linear,
        "pjc": pjc,
    }


def operator_analysis(derivation: Mapping | None = None) -> dict:
    derivation = derivation or derive_jump_classes()
    classes = derivation["classes"]
    jumps = {
        name: [event["delta"][index] for event in classes]
        for index, name in enumerate(OBSERVABLES)
    }
    span = discover_diagonal_commutator_span(
        event_classes=[event["id"] for event in classes], jumps=jumps
    )
    operator = emgf_operator_record(
        observable_names=OBSERVABLES,
        event_classes=classes,
        variables=VARIABLES,
        jump_projection=REMY_JUMP_PROJECTION,
    )
    # The monomials are generated from the jump vectors; this is also the
    # direct consistency check for the variable substitution used in Eq. 58.
    monomials = []
    variable_names = ("x_E", "x_P1", "x_P2", "x_P3")
    for event in classes:
        numerator = []
        denominator = []
        for variable, power in zip(variable_names, event["delta"]):
            target = numerator if power > 0 else denominator
            magnitude = abs(power)
            if magnitude == 1:
                target.append(variable)
            elif magnitude:
                target.append(f"{variable}^{magnitude}")
        factor = "*".join(numerator) or "1"
        if denominator:
            factor += "/" + "*".join(denominator)
        monomials.append({"event_class": event["id"], "factor": factor})
    return {
        "span": span.record(),
        "emgf_operator": operator,
        "variable_substitution_monomials": monomials,
        "source_comparison": {
            "eq57g": "generated coefficients agree when collected by derivatives",
            "eq58": (
                "the printed monomials do not follow literally from omega_i=ln(x_i); "
                "the generated factors above follow from the derived jumps"
            ),
        },
    }
