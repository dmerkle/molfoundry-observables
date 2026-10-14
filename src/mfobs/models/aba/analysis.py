"""Run path support discovery, event extraction and linear observable closure."""
from __future__ import annotations

import mfobs.models.aba.events as discovery_events
from mfobs.outputs import support_record as export_support
from hashlib import sha256
import json
from time import perf_counter
from numbers import Integral
from mfobs.core.analysis_contract import (
    AnalysisContract, NormalizedRule, ObservableDeclaration, StateLanguageContract,
)
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.models.aba import counting as discovery


def _contract(model, support, sources, selected_queries):
    gated = model["gate"]
    normalized = tuple(NormalizedRule(
        source.identifier, source.rate, "conditional-DPO",
        source.rule.getGMLString(),
        ("NoEdge(GateSource,GateTarget)",) if gated and source.kind in ("append", "bind") else (),
        boundary=source.kind == "feed",
    ) for source in sources)
    meanings = {
        "motif": "occurrences in the oriented chain backbone",
        "suffix": "indicator that the unique chain ends in this word",
        "pool": "number of free monomer copies",
        "gate": "indicator that the unique NoEdge gate is open",
        "total": "free plus incorporated copies of the selected symbol",
    }
    observables = tuple(ObservableDeclaration(
        name, meanings[pattern.kind],
        "selected" if name in selected_queries else "auxiliary",
        ("user-selected query" if name in selected_queries else
         "discovered graph factor; exact origin and composition trace in support_discovery"),
        pattern.rule.getGMLString(),
    ) for name, pattern in support.patterns.items())
    language = StateLanguageContract(
        f"terminal-word-{model['identifier']}",
        "one rooted oriented chain, arbitrary finite A/B monomer pools, persistent feeds"
        + (", and one Boolean NoEdge gate" if gated else ""),
        True,
        ("exactly one rooted directed path and one End marker",
         "typed projection uses one End and saturated two-port Next/Terminal incidence nodes",
         "free monomer and incorporated-repeat matches are disjoint",
         "symbols retain their oriented backbone order",
         * (("unique GateSource/GateTarget pair; X=1 iff its bound edge is absent",)
            if gated else ())),
        ("reject nonfunctional, cyclic and multiply rooted projected overlaps",
         "erase a DPO deletion interface only for saturated deleted incidence nodes",
         "factor only components with disjoint match images in this language",
         "refuse unsupported constraints and graph shapes"),
        "Feed adds a free symbol-tagged position; insertion attaches it to the unique "
        "path end; gate rules toggle only the GateSource--GateTarget relation.",
        "word, free-pool counts" + (", and Boolean gate state" if gated else ""),
    )
    return AnalysisContract(
        f"aba-{model['identifier']}", model["title"], normalized, language,
        observables, tuple(selected_queries), tuple(selected_queries),
        {source.identifier: "1 per injective match" for source in sources},
    )


def analyze(inputs, *, linear_limit=24,
                               diagnostic_timings=False,
                               support_pattern_limit=64,
                               support_expansion_limit=64,
                               record_compositions=False,
                               run_linear=True):
    started = perf_counter()
    timings = {}

    def mark(name, beginning):
        if diagnostic_timings:
            timings[name] = perf_counter() - beginning

    state_space = inputs.omega
    alphabet, gated = state_space.definition["alphabet"], state_space.definition["gated"]
    queries = {}
    motifs = []
    for name, observable in inputs.observables.items():
        terms = observable.identity_terms()
        if len(terms) != 1 or terms[0][0] != 1:
            raise ValueError("ABA requires one unit-coefficient motif per query")
        reduced = discovery.reduce_input(terms[0][1])
        if len(reduced.factors) != 1:
            raise ValueError("ABA query must be a single path motif")
        pattern = next(iter(reduced.factors.values()))
        if pattern.name != name or pattern.kind != "motif":
            raise ValueError("ABA query must use its canonical motif name")
        queries[name] = pattern
        motifs.append(pattern.metadata["word"])
    query_names = tuple(queries)
    model = {"identifier": "-".join(motifs) + ("-gated" if gated else ""),
             "title": "Growing paths", "alphabet": list(alphabet),
             "motifs": motifs, "gate": gated}
    effective = {"schema_version": 1, "kind": "path-analysis-configuration",
                 "model": model, "linear_dimension_limit": linear_limit}
    if not query_names:
        raise ValueError("support discovery requires at least one selected query")
    if run_linear and linear_limit < len(query_names):
        raise ValueError(f"linear dimension limit {linear_limit} is smaller than the "
                         f"{len(query_names)} selected linear-discovery seeds")
    for value, minimum, label in (
        (support_pattern_limit, len(query_names), "support pattern"),
        (support_expansion_limit, 0, "support expansion"),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
            raise ValueError(f"{label} limit must be an integer at least {minimum}")
    effective.update({
        "representation": "directed paths", "selected_queries": list(query_names),
        "support_pattern_limit": support_pattern_limit,
        "support_expansion_limit": support_expansion_limit,
    })
    effective["sha256"] = sha256(json.dumps(
        effective, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()

    phase = perf_counter()
    linear = None
    if run_linear:
        rules = tuple(inputs.rules.values())
        linear = discover_linear_observables(rules, queries,
            counter=discovery.PathCounter(state_space=state_space), maximum_dimension=linear_limit)
        mark("linear_discovery", phase)
    phase = perf_counter()
    support, sources = discovery.derive(
        motifs, alphabet=alphabet, gated=gated,
        maximum_patterns=support_pattern_limit,
        maximum_expansions=support_expansion_limit,
        record_compositions=record_compositions, state_space=state_space,
    )
    mark("support_discovery", phase)
    phase = perf_counter()
    support_record = export_support(support, sources, events=discovery_events.compile_uniform_events(support, sources))
    event_record = support_record["events"]
    contract = _contract(model, support, sources, query_names)
    definition_sha = contract.digest
    projection = {
        "kind": "typed terminal-path rules",
        "structural_class_count": len(event_record["events"]),
        "source_rule_count": len(sources),
        "source_origin": "query-independent depth-zero typed constructors",
        "molecular_correspondence": "separate atom-resolved analysis and comparison",

    }
    composition = {
        "status": "PASS" if event_record["status"] == "FINITE_EVENTS" else event_record["status"],
        "backend": "MolFoundry DPO common-overlap composition with discovered graph factors",
        "jump_projection": "DPO",
        "observable_coordinates": list(support.patterns),
        "event_classes": event_record["events"],
        "projection": projection,
        "derivation_trace": "support_discovery.trace",
    }
    result = {
        "schema_version": 1, "variant": model["identifier"],
        "configuration": effective, "contract": contract.record(),
        "support_discovery": support_record, "composition": composition,
        "linear": linear,
        "assumptions_evidence": {
            "schema_version": 1,
            "input_definition_sha256": definition_sha,
            "configuration_sha256": effective["sha256"],
            "scope": "supported terminal-path constructors on all finite words and pools",
            "state_language": contract.state_language.record(),
            "support_status": support.status,
            "event_status": event_record["status"],
            "uniform_match_argument": event_record.get("uniform_match_obligations", []),
            "composition_trace": "support_discovery.trace",

            "finite_enumeration_used_by_analysis": False,
        },
    }
    if linear is not None:
        linear["input_definition_sha256"] = definition_sha
    if event_record["status"] != "FINITE_EVENTS":
        result["status"] = (support.status if support.status.startswith("UNKNOWN")
                            else event_record["status"])
        reason = event_record.get("reason") or support.reason or event_record["status"]
        result["pjc"] = {
                "schema_version": 1, "status": "UNKNOWN_PREREQUISITE",
                "input_definition_sha256": definition_sha, "reason": reason,
                "prerequisite_status": result["status"],
        }
    else:
        model = discovery_events.polynomial_jump_model(
            support, sources, name=model["identifier"], input_definition_sha256=definition_sha,
            extracted=event_record,
        )
        # The usual composition record exposes the quotient-reduced events;
        # raw discovery activities and predicates remain in support_discovery.
        composition["event_classes"] = [{
            "identifier": event.identifier, "source_rule": event.rule,
            "rate": event.channel, "predicate": event.predicate,
            "jump": {n: value for n, value in zip(model.ring.variables, event.delta) if value},
            "activity": model.ring.record(event.propensity),
        } for event in model.events]
        composition["jump_model_provenance"] = model.provenance.record()
        mark("event_compilation", phase)
        phase = perf_counter()
        result["pjc"] = analyze_polynomial_jump_closure(model)
        mark("polynomial_jump_closure", phase)
        result["status"] = "ANALYZED"
    if diagnostic_timings:
        timings["total"] = perf_counter() - started
        result["diagnostics"] = {
            "kind": "wall-clock phase timings",
            "scope": "diagnostic only; support and linear worklist budgets are reported separately",
            "seconds": timings,
        }
    return result


def word_coordinates(record: dict, word: str, pools: dict[str, int], opened=True) -> dict[str, int]:
    """Evaluate the *discovered* pattern meanings on one user-selected state."""
    patterns = record.get("analysis", record)["support_discovery"]["patterns"]
    if not isinstance(word, str) or not isinstance(opened, bool):
        raise ValueError("Supply a string word and a Boolean gate state")
    alphabet = {p.get("word") for p in patterns.values() if p["kind"] == "pool"}
    if set(pools) != alphabet or set(word) - alphabet:
        raise ValueError("Word symbols and pool keys must match the recorded free-pool alphabet")
    if any(not isinstance(v, Integral) or isinstance(v, bool) or v < 0 for v in pools.values()):
        raise ValueError("Pool counts must be nonnegative integers")
    values = {}
    for name, pattern in patterns.items():
        motif = pattern.get("word", "")
        kind = pattern["kind"]
        if kind == "pool":
            value = pools[motif]
        elif kind == "total":
            value = pools[motif] + word.count(motif)
        elif kind == "motif":
            if not motif:
                raise ValueError("Empty motif")
            value = sum(word[i:i + len(motif)] == motif for i in range(len(word) - len(motif) + 1))
        elif kind == "suffix":
            value = int(word.endswith(motif))
        elif kind == "gate":
            value = int(opened)
        else:
            raise ValueError(f"Unknown recorded pattern kind {kind!r}")
        values[name] = int(value)
    if set(values) != set(record.get("analysis", record)["pjc"]["observable_coordinates"]):
        raise ValueError("Discovered pattern meanings and event coordinates disagree")
    return values
