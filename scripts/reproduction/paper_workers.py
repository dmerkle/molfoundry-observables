"""Adapters for the paper claims; no supplied target dynamics.

Each function invokes an existing analysis or validation entry point. Expected
numbers live only in the manifest and are checked after analysis has completed.
"""
from mfobs.applications import MolecularROP
from pathlib import Path


def aldol(stage, spec):
    from scripts.reproduction.reproduce_paper import input_identity, selected_manifest
    from examples.tcs.compare_conditions import compare_conditions, HERE
    report = compare_conditions(Path(spec["_output_dir"]) / "compositions")
    observed = {}
    for name, data in report["variants"].items():
        observed[name] = {
            "after": data["orders"]["reaction_then_query"]["filtered"],
            "before": data["orders"]["query_then_reaction"]["filtered"],
            "signed_terms": data["nonzero_signed_terms"],
        }
    observed["removed_after"] = sorted(
        row["plain_class"] for row in report["comparison"]["removed_contributions"]["reaction_then_query"])
    observed["removed_before"] = len(report["comparison"]["removed_contributions"]["query_then_reaction"])
    observed["r180_coefficients"] = {
        name: sorted(row["net_coefficient"] for row in rows)
        for name, rows in report["comparison"]["selected_r180"].items()
    }
    return {"result": report, "observed": observed,
            "inputs": input_identity([selected_manifest(spec), *sorted(HERE.glob("*.gml")),
                                      HERE/"chemistry.json", HERE/"article-aliases.json"])}


def path_query(stage, spec):
    from scripts.reproduction.reproduce_paper import analysis_summary, input_identity, selected_manifest, PROJECT, json_digest
    from examples.aba.run import run as run_aba
    query = {"motif": spec["motif"], "alphabet": ["A", "B"], "gate": spec.get("gated", False),
             "linear_limit": 24, "support_pattern_limit": 64, "support_expansion_limit": 64}
    analysis = run_aba(query, diagnostic_timings=True)
    from mfobs.checks.terminal_paths import validate_result
    finite = validate_result(analysis)
    row = analysis_summary(stage, analysis, finite)
    support = analysis["support_discovery"]
    row.update(support_status=support["status"], coordinate_count=len(support["patterns"]),
               suffixes=sorted(n for n, p in support["patterns"].items() if p["kind"] == "suffix"),
               source_rules_sha256=json_digest(support["source_rules"]))
    return {"analysis": analysis, "analysis_row": row, "observed": row,
            "query": query, "finite_validation": finite,
            "inputs": input_identity([selected_manifest(spec), PROJECT/"examples/aba/run.py"])}


def molecular(stage, spec):
    from scripts.reproduction.reproduce_paper import analysis_summary, input_identity, selected_manifest, PROJECT
    from mfobs.models.molecular_rop import inputs as rop
    selection = spec.get("preset") or rop.configured_variant(motifs=("ABA",), gated=spec.get("gated", False))
    from mfobs import analyze
    backend = MolecularROP(selection)
    analysis = analyze(backend.specification(), backend=backend, linear_limit=24, molecular_overlap_limit=spec.get("overlap_limit", 2000000), diagnostic_timings=True).record()
    from mfobs.checks.molecular_rop import validate_result
    checks = validate_result(analysis)
    finite = checks["molecular_graph_validation"]
    row = analysis_summary(stage, analysis, finite)
    row.update(status=analysis["status"], coordinate_count=len(analysis["support_discovery"]["patterns"]),
               word_check=checks["word_oracle_comparison"]["status"],
               typed_check=checks["typed_reference_comparison"]["status"])
    if "overlap_limit" in spec:
        support = analysis["support_discovery"]
        failures = [r.get("minimal_reproducer", {}) for r in support["trace"]]
        row.update(pending_preserved=bool(support["pending"]),
                   operands_preserved=any(r.get("query_gml") and r.get("source_gml") for r in failures),
                   operator_absent="operator" not in analysis["pjc"])
    return {"analysis": analysis, "analysis_row": row, "observed": row,
            "finite_validation": finite,
            "inputs": input_identity([selected_manifest(spec), PROJECT/"src/mfobs/models/molecular_rop/inputs.py",
                                      PROJECT/"src/mfobs/models/molecular_rop/counting.py"])}


def eight_slots(stage, spec):
    from scripts.reproduction.reproduce_paper import input_identity, selected_manifest, PROJECT
    from mfobs import EightSlots, analyze
    from scripts.check_eight_slots import verify
    output = Path(spec["_output_dir"]) / "symbolic"
    backend = EightSlots()
    result = analyze(backend.specification(), backend=backend,
                     max_expansions=spec["max_expansions"]).record()
    # No oracle data enters analysis. Only its computed result is passed to
    # the independent checker, which then enumerates the three topologies.
    checked = verify(output, result=result)
    observed = {"closed": result["closed"], "coordinate_count": len(result["basis"]),
                "pending": result["pending"], "unknown": result["unknown"],
                "synthesis_state_evaluations": result["synthesis_state_evaluations"],
                "synthesis_orbit_evaluations": result["synthesis_orbit_evaluations"],
                "verification_passed": checked["passed"], "rank": checked["verified_count_function_rank"],
                "states": sum(r["occupancy_assignments_checked"] for r in checked["families"].values()),
                "pattern_checks": checked["independent_raw_pattern_vs_polynomial_checks"],
                "generator_rows": checked["exact_rate_separated_generator_rows"]}
    return {"result": result, "finite_validation": checked, "observed": observed,
            "inputs": input_identity([selected_manifest(spec), PROJECT/"examples/eight-slots/model.json"])}


def dna_walker(stage, spec):
    import json
    from mfobs import DNAWalker, analyze
    from mfobs.checks.dna_walker import validate, verified_exports, compare_incidence
    from scripts.reproduction.reproduce_paper import input_identity, selected_manifest, PROJECT
    input_path = PROJECT / "examples/dna-walker/input.json"
    inputs = json.loads(input_path.read_text(encoding="utf-8"))
    backend = DNAWalker()
    result = analyze(backend.specification(), backend=backend).record()
    evidence = PROJECT / "notebooks/visual/evidence/dna-walker"
    archived = json.loads((evidence / "result.json").read_text(encoding="utf-8"))
    provenance = json.loads((evidence / "provenance.json").read_text(encoding="utf-8"))
    original, incidence = verified_exports(evidence / "recovered", archived, provenance)
    assert inputs == {"rules": original["rules"], "query": original["extended_query"]}
    comparison = compare_incidence(result, incidence)
    finite = validate(result, inputs["rules"], {
        "E": original["extended_query"], "Craw": original["compressed_query"]})
    observed = {"status": result["status"], "dimension": len(result["basis"]),
                "fragger_contributions_compared": len(comparison),
                "fragger_mismatches": sum(not r["Equal"] for r in comparison),
                "validation_status": finite["status"],
                "hosts": finite["hosts"],
                "generator_images_checked": finite["generator_images_checked"],
                "synthesis_state_evaluations": result["synthesis_state_evaluations"]}
    return {"result": result, "observed": observed, "finite_validation": finite,
            "fragger_comparison": comparison,
            "comparison_scope": "Fresh MolFoundry discovery and literal host checks; retained FRAGGER execution outputs.",
            "inputs": input_identity([selected_manifest(spec), input_path,
                evidence / "result.json", evidence / "provenance.json",
                *sorted((evidence / "recovered").glob("*.json"))])}


def context_nac(stage, spec):
    from mfobs import Chlorination, analyze
    from scripts.reproduction.reproduce_paper import input_identity, selected_manifest, PROJECT
    backend = Chlorination(spec["variant"])
    analysis = analyze(backend.specification(), backend=backend).record()
    from mfobs.checks.chlorination import validate_all_masks
    finite = validate_all_masks(backend, analysis)
    row = {"stage": stage, "status": analysis["status"],
           "pjc_status": analysis["pjc"]["status"],
           "event_count": len(analysis["events"]),
           "coordinate_count": len(analysis["pjc"]["coordinates"]),
           "finite_status": finite["status"],
           "joint_rate_checks": finite["molfoundry_joint_rate_checks"],
           "growth_rate_checks": finite["molfoundry_growth_rate_checks"]}
    return {"analysis": analysis, "analysis_row": row, "observed": row,
            "finite_validation": finite,
            "inputs": input_identity([selected_manifest(spec),
                PROJECT/"src/mfobs/models/chlorination/inputs.py"])}


WORKERS = {"context-nac": context_nac, "aldol-conditions": aldol, "path-query": path_query,
           "molecular-rop": molecular, "eight-slots": eight_slots, "dna-walker": dna_walker}
