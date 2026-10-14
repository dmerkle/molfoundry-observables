"""Uniform-match events on growing paths using discovered suffix and gate predicates."""
from __future__ import annotations
from itertools import product
import sympy as sp
from mfobs.core.pattern_support import Support, UnsupportedReduction
from mfobs.models.aba.inputs import source_rules


def consistent_suffix_assignment(assignment, patterns) -> bool:
    # All positive suffixes must be suffixes of one longest positive suffix.
    # That word itself witnesses consistency unless a negative suffix forbids it.
    positives = [patterns[n].metadata["word"] for n, bit in assignment.items()
                 if bit and patterns[n].kind == "suffix"]
    witness = max(positives, key=len, default="")
    return all(witness.endswith(patterns[n].metadata["word"]) == bool(bit)
               for n, bit in assignment.items() if patterns[n].kind == "suffix")


def compile_uniform_events(support: Support, sources) -> dict:
    """Finite Boolean partition of composition-derived state-dependent jumps.

    Additional all-language obligation: at a fixed valid host, all matches of
    one source rule have isomorphic successors. Unique terminal/gate/feed and
    interchangeable free tokens establish it for these constructors. The caller
    cannot infer it from polynomial support closure. This method accepts only
    the literal source constructors below, checking conditional canonical keys.
    """
    if support.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status": "UNKNOWN_SUPPORT", "events": []}
    alphabet = tuple(r.symbol for r in sources if r.kind == "feed")
    try:
        expected = source_rules(alphabet, any(r.kind == "bind" for r in sources))
    except ValueError as exc:
        return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "reason": str(exc), "events": []}
    if len(sources) != len(expected) or any(
        (r.identifier, r.kind, r.symbol) != (e.identifier, e.kind, e.symbol)
        or r.rule.canonKey != e.rule.canonKey
        for r, e in zip(sources, expected)
    ):
        return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "events": []}
    return partition_uniform_events(support, sources,
        scope="all finite words and finite pools in the declared terminal language",
        uniform_match_argument="unique terminal/feed/gate and isomorphic free-token choices")


def partition_uniform_events(support, sources, *, scope, uniform_match_argument):
    """Terminal-language partition; the event algebra lives in the core."""
    from mfobs.core.uniform_events import EventClass, extract_uniform_events

    def terminal_partition(source, jumps):
        dependencies = sorted({str(s) for value in jumps.values() for s in value.free_symbols})
        if any(support.patterns[n].kind not in ("suffix", "gate") for n in dependencies):
            raise UnsupportedReduction("unbounded state-dependent jumps")
        for values in product((0, 1), repeat=len(dependencies)):
            assignment = dict(zip(dependencies, values))
            if consistent_suffix_assignment(assignment, support.patterns):
                indicator = sp.prod(sp.Symbol(n) if bit else 1-sp.Symbol(n)
                                    for n, bit in assignment.items())
                yield EventClass({sp.Symbol(n): bit for n, bit in assignment.items()},
                                 indicator, assignment)

    return extract_uniform_events(support, sources, partition=terminal_partition,
        scope=scope, uniform_match_argument=uniform_match_argument,
        partition_argument="Boolean predicates on DPO-discovered support patterns")


def polynomial_jump_model(support: Support, sources, *, name: str,
                          input_definition_sha256: str = "", extracted=None, event_compiler=None,
                          backend="query-driven terminal-path support discovery",
                          jump_projection="DPO with checked saturated deletion interfaces"):
    """Adapt certified events to the existing exact linear/EMGF consumers.

    Boolean products are language identities, not statistical independence:
    two compatible suffix indicators multiply to the longer suffix, while
    incompatible suffixes cannot both hold. A gate is independently Boolean.
    This constructor refuses incomplete support or an unproved event route.
    """
    from mfobs.core.polynomial_jump_closure import (
        JumpModelProvenance, PolynomialRing,
    )
    if extracted is None:
        extracted = (event_compiler or compile_uniform_events)(support, sources)
    if extracted["status"] != "FINITE_EVENTS":
        raise UnsupportedReduction(extracted.get("reason", extracted["status"]))
    names = tuple(support.patterns)
    suffixes = {n: p.metadata["word"] for n, p in support.patterns.items()
                if p.kind == "suffix"}
    products = {}
    for index, left in enumerate(suffixes):
        for right in tuple(suffixes)[index + 1:]:
            lw, rw = suffixes[left], suffixes[right]
            products[left, right] = (left if lw.endswith(rw) else
                                     right if rw.endswith(lw) else None)
    ring = PolynomialRing(names, boolean_variables=[
        n for n, p in support.patterns.items() if p.kind in ("suffix", "gate")
    ], boolean_products=products)

    certificate = {
        "status": "VERIFIED_ALL_STATES",
        "scope": extracted["scope"],
        "argument_status": "DECLARED_CONSTRUCTOR_ARGUMENT_WITH_CHECKED_SOURCE_IDENTITIES",
        "partition": extracted["event_partition_origin"],
        "source_identity_check": "conditional canonical keys, identifiers, kinds and symbols",
        "uniform_match_obligations": extracted["uniform_match_obligations"],
        "nonnegative_activities": "raw injective source-match counts times Boolean class indicators",
        "state_preservation": "each supported source feeds a free token, appends it, or toggles the unique gate",
        "proof_scope": extracted.get("proof_scope", "supported literal constructors on the declared terminal-path language; not inferred from finite validation"),
    }
    from mfobs.core.uniform_events import compile_event_model
    return compile_event_model(
        extracted, name=name, ring=ring, certificate=certificate,
        input_definition_sha256=input_definition_sha256,
        provenance=JumpModelProvenance(
            backend=backend,
            event_partition_origin=extracted["event_partition_origin"],
            jump_origin="unrefined DPO generator images divided by source activity under the uniform-match argument",
            activity_origin="exact projected input counts times discovered Boolean predicates",
            jump_projection=jump_projection,
            composition_scope="selected and recursively discovered patterns against unrefined source rules",
            all_jumps_composition_derived=True, supplied_delta_table=False,
        ),
    ).model
