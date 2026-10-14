"""Run atom-resolved support discovery, event extraction and linear closure."""
from __future__ import annotations

import mfobs.models.molecular_rop.events as discovery_events
import mfobs.models.molecular_rop.inputs as discovery_inputs
from mfobs.outputs import support_record as export_support
from mfobs.provenance import source_hashes
from hashlib import sha256
from importlib.metadata import version
import json
from time import perf_counter
from mfobs.core.analysis_contract import AnalysisContract, NormalizedRule, ObservableDeclaration, StateLanguageContract
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.molecular_rop import counting as discovery


def _contract(selected,support,sources,query_names):
    language = StateLanguageContract(
        f"direct-molecular-rop-{selected.identifier}",
        "one benzyl-initiated unbranched lactone chain, arbitrary finite free A/B lactone pools, persistent feed markers"
        + (" and one asymmetric Pd-N gate" if selected.gated else ""),
        True,
        ("constructor-defined molecular graphs with explicit hydrogens",
         "one terminal OH; oriented complete repeat and free-lactone matches have unit normalization",
         "free pools, chain, gate and feed markers are disjoint host components"),
        ("compose chemical rules and chemical queries before any sequence interpretation",
         "retain all common overlaps, including disconnected and empty overlaps",
         "reject only certified valence/bond violations and patterns absent from the complete local-context cover",
         "reduce only represented whole-repeat/pool/end patterns and molecular Pd-N NoEdge",
         "do not erase DPO vertex-deletion interfaces; return UNKNOWN"),
        "Feed adds a complete free lactone; atom-preserving alcoholysis extends the unique terminal OH; gate rules toggle only Pd-N coordination. These rules preserve the constructor family.",
        "derived pattern counts; sequence meanings are assigned only during post-composition reduction",
        finite_enumeration_role="finite local-context cover certifies zero connected patterns for all chains",
    )
    return AnalysisContract(
        f"molecular-rop-{selected.identifier}",selected.title,
        tuple(NormalizedRule(r.identifier,r.rate,"conditional-DPO",r.rule.getGMLString(),boundary=r.kind=="feed") for r in sources),
        language,
        tuple(ObservableDeclaration(n, p.kind+" molecular pattern count", "selected" if n in query_names else "auxiliary",
            "selected chemical query" if n in query_names else "molecular component discovered in a reduced chemical composition",
            p.rule.getGMLString()) for n,p in support.patterns.items()),
        tuple(query_names),tuple(query_names),{r.identifier:"1 per injective atom map" for r in sources},
    )


def analyze(selection,*,linear_limit=24,diagnostic_timings=False,support_pattern_limit=64,support_expansion_limit=64,molecular_overlap_limit=2000000,state_space=None,run_linear=True):
    started=perf_counter(); selected=rop.configuration(selection)
    motifs=tuple("".join(w) for w in selected.motifs)
    if selected in rop.VARIANTS.values(): motifs=selected.alphabet+motifs
    query_names=tuple(("N_" if len(w)==1 else "C_")+w for w in motifs)
    if not query_names: raise ValueError("at least one molecular query is required")
    for value,minimum,label in [(linear_limit,len(query_names) if run_linear else 1,"linear"),(support_pattern_limit,len(query_names),"support pattern"),(support_expansion_limit,0,"support expansion"),(molecular_overlap_limit,0,"molecular overlap")]:
        if not isinstance(value,int) or isinstance(value,bool) or value<minimum:
            raise ValueError(f"{label} budget must be an integer at least {minimum}")
    effective=rop.effective_configuration(selected,linear_limit=linear_limit)
    effective.pop("sha256")
    effective.update(support_method="molecular",selected_queries=list(query_names),support_pattern_limit=support_pattern_limit,
        support_expansion_limit=support_expansion_limit,maximum_local_context_depth=8,maximum_candidate_overlaps_per_order=molecular_overlap_limit)
    effective["sha256"]=sha256(json.dumps(effective,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    linear = None
    if run_linear:
        sources = discovery.source_rules(selected, state_space=state_space)
        queries = {q.name: q for q in map(discovery.query_pattern, motifs)}
        linear = discover_linear_observables(sources, queries,
            counter=discovery.MolecularReducer(selected, maximum_candidates=molecular_overlap_limit, state_space=state_space),
            maximum_dimension=linear_limit)
    linear_seconds = perf_counter()-started
    support_started=perf_counter()
    support,sources=discovery.derive(selected,motifs,maximum_patterns=support_pattern_limit,maximum_expansions=support_expansion_limit,maximum_candidates=molecular_overlap_limit,state_space=state_space)
    support_seconds=perf_counter()-support_started
    record=export_support(support,sources,events=discovery_events.compile_uniform_events(support,sources)); events=record["events"]
    contract=_contract(selected,support,sources,query_names)
    result={"schema_version":1,"variant":selected.identifier,"configuration":effective,"contract":contract.record(),
        "support_discovery":record,
        "linear":linear,
        "composition":{"status":"PASS" if events["status"]=="FINITE_EVENTS" else events["status"],
            "backend":"direct atom-resolved molecular DPO composition","jump_projection":"DPO",
            "observable_coordinates":list(support.patterns),"event_classes":events["events"],
            "projection":{"kind":"no typed projection before composition","source_rule_count":len(sources),
                "structural_class_count":len(events["events"]),"source_origin":"literal atom-resolved molecular_rop.inputs.rules constructors"},
            "derivation_trace":"support_discovery.trace"},
        "assumptions_evidence":{"schema_version":1,"input_definition_sha256":contract.digest,
            "configuration_sha256":effective["sha256"],"state_language":contract.state_language.record(),
            "scope":"constructor-defined molecular family only","support_status":support.status,"event_status":events["status"],
            "uniform_match_argument":events.get("uniform_match_obligations",[]),
            "local_context_cover_used_for_zero_certificates":True,"finite_state_interpolation_used":False,
            "composition_restrictions":discovery_inputs.forbidden_specification(selected,maximum_candidates=molecular_overlap_limit)},
        "implementation_provenance":{"molfoundry_version":version("molfoundry"),"sympy_version":version("sympy"),"networkx_version":version("networkx"),
            "source_sha256":source_hashes()},
    }
    if linear is not None: linear["input_definition_sha256"]=contract.digest
    if events["status"]!="FINITE_EVENTS":
        result["status"]=support.status if support.status.startswith("UNKNOWN") else events["status"]
        reason=events.get("reason") or support.reason or events["status"]
        result["pjc"]={"status":"UNKNOWN_PREREQUISITE","reason":reason,"prerequisite_status":result["status"],"input_definition_sha256":contract.digest}
    else:
        model=discovery_events.polynomial_jump_model(support,sources,name="molecular-"+selected.identifier,input_definition_sha256=contract.digest,extracted=events)
        result["composition"]["event_classes"]=[{"identifier":e.identifier,"source_rule":e.rule,"rate":e.channel,"predicate":e.predicate,
            "jump":{n:v for n,v in zip(model.ring.variables,e.delta) if v},"activity":model.ring.record(e.propensity)} for e in model.events]
        result["composition"]["jump_model_provenance"]=model.provenance.record()
        result["pjc"]=analyze_polynomial_jump_closure(model)
        result["status"]="ANALYZED"
    if diagnostic_timings:
        result["diagnostics"]={"kind":"wall-clock phase timings","scope":"fresh direct chemical composition only; not historical typed-path timings",
            "seconds":{"linear_discovery":linear_seconds,"molecular_support_discovery":support_seconds,"total":perf_counter()-started}}
    return result
