"""Uniform molecular source matches and the discovered suffix event partition."""
from __future__ import annotations

import mfobs.models.aba.events as terminal_events
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.molecular_rop.inputs import source_rules


def compile_uniform_events(support, sources):
    """Check chemical constructor identities before using shared partition algebra."""
    if support.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status": "UNKNOWN_SUPPORT", "events": []}
    alphabet = tuple(r.symbol for r in sources if r.kind == "feed")
    gated = any(r.kind == "bind" for r in sources)
    try:
        selection = rop.RopConfiguration("molecular-source-check", "Molecular source check", alphabet, gated, (), "source identity check")
        expected = source_rules(selection)
    except (ValueError, KeyError) as exc:
        return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "reason": str(exc), "events": []}
    if len(sources) != len(expected) or any(
        (a.identifier,a.kind,a.symbol,a.rate,a.rule.canonKey) !=
        (b.identifier,b.kind,b.symbol,b.rate,b.rule.canonKey)
        for a,b in zip(sources,expected)
    ):
        return {"status": "UNKNOWN_UNIFORM_MATCH_OBLIGATION", "reason": "chemical rules differ from the justified constructor family", "events": []}
    result = terminal_events.partition_uniform_events(support,sources,
        scope="all constructor-defined finite molecular chains and free-lactone pools, with optional unique Pd-N gate",
        uniform_match_argument="unique molecular terminal OH and oriented free-lactone match; different monomer copies yield isomorphic successors; unique feed markers and asymmetric Pd-N gate")
    result["event_partition_origin"] = "Boolean predicates on molecular-composition-discovered pattern counts"
    result["proof_scope"] = "literal molecular constructors, induced whole-unit pattern identities, disjoint component factorization and complete local-context zero certificates; bounded host validation is separate"
    return result


def polynomial_jump_model(support,sources,*,name,input_definition_sha256="",extracted=None):
    return terminal_events.polynomial_jump_model(support,sources,name=name,input_definition_sha256=input_definition_sha256,extracted=extracted,
        event_compiler=compile_uniform_events, backend="direct molecular ROP composition and pattern discovery",
        jump_projection="DPO retaining the interface; only no-vertex-deletion projected inputs reduced")
