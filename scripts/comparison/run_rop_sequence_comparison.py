#!/usr/bin/env python3
"""Bounded external FRAGGER diagnostic for the positive ungated ROP word graph.

The runner exports the eight concrete directed graph rules actually supplied
to FRAGGER, validates their declared word semantics on a finite nonvacuous
domain, and records two fresh bounded worklists. Every graph factor supported
by the conservative language mapping is compared with the current
composition-derived MolFoundry generator. A process failure or worklist cap is
never reported as a nonclosure claim.
"""
from __future__ import annotations

import mfobs.checks.rop_projection as projection_checks

import argparse
import hashlib
import json
from collections import Counter
from itertools import product
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mfobs.core import engine
from mfobs.core.source_identity import (
    FRAGGER_SNAPSHOT,
    require_source_identity,
    source_identity,
)
from mfobs.core.polynomial_jump_closure import generator_action

COMMIT = "0d9595c5d6e8918fc605f7d65778568a74f5e2d4"
JAR_HASHES = {
    "scala-compiler-2.13.2.jar": "b04da83111decff5b36e87c74d73f5e2231b63ce605f2602e862bae917d7761a",
    "scala-library-2.13.2.jar": "6a1c9846feb8f5d7d33b91e61e056d272215a79747a2b7b5e5ece690bf151d5d",
    "scala-reflect-2.13.2.jar": "21a7288874f186d6049e8b14a986011ffd96cb587f206eccf946624f2d476815",
}
QUERIES = ("N_A", "N_B", "C_AA", "C_AB", "C_BA", "C_BB", "C_ABA")


def command(args, cwd: Path, timeout: int, *, check: bool = True):
    return subprocess.run(
        [str(value) for value in args],
        cwd=cwd,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=check,
    )


def external_worklist(
    args, cwd: Path, timeout: int, log_path: Path, budget: int
) -> dict:
    """Run one bounded process and preserve timeout, execution, and parse failures."""
    started = time.perf_counter()
    base = {
        "budget": budget,
        "command": [str(value) for value in args],
        "timeout_seconds": timeout,
        "log": str(log_path),
    }
    try:
        run = command(args, cwd, timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        stdout = (
            exc.stdout.decode()
            if isinstance(exc.stdout, bytes)
            else (exc.stdout or "")
        )
        stderr = (
            exc.stderr.decode()
            if isinstance(exc.stderr, bytes)
            else (exc.stderr or "")
        )
        log_path.write_text(stdout + stderr, encoding="utf-8")
        return {
            **base,
            "process_seconds": time.perf_counter() - started,
            "returncode": None,
            "status": "TIMEOUT",
            "json": None,
            "error": f"process exceeded {timeout} seconds",
        }
    log_path.write_text(run.stdout + run.stderr, encoding="utf-8")
    result = {
        **base,
        "process_seconds": time.perf_counter() - started,
        "returncode": run.returncode,
        "json": None,
    }
    if run.returncode:
        return {
            **result,
            "status": "EXECUTION_ERROR",
            "error": f"process exited with status {run.returncode}",
        }
    # Upstream creates fresh node IDs from arbitrary Unicode characters.  U+2028
    # and similar code points are valid inside JSON strings but `splitlines()`
    # treats them as record boundaries, so only ASCII LF delimits harness output.
    json_lines = [
        line
        for line in run.stdout.split("\n")
        if line.startswith("FRAGGER_JSON=")
    ]
    if len(json_lines) != 1:
        return {
            **result,
            "status": "PARSE_ERROR",
            "error": f"expected one FRAGGER_JSON record, found {len(json_lines)}",
        }
    try:
        payload = json.loads(json_lines[0].split("=", 1)[1])
    except json.JSONDecodeError as exc:
        return {**result, "status": "PARSE_ERROR", "error": str(exc)}
    return {**result, "status": "PASS", "json": payload}


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def path_graph(word: tuple[str, ...], pools: tuple[int, int] = (0, 0)) -> dict:
    """Complete declared host serialization, including isolated free vertices."""
    nodes = [{"id": "r", "label": "R"}]
    nodes += [{"id": f"w{i}", "label": label} for i, label in enumerate(word)]
    nodes += [{"id": "t", "label": "T"}]
    nodes += [{"id": f"fa{i}", "label": "F_A"} for i in range(pools[0])]
    nodes += [{"id": f"fb{i}", "label": "F_B"} for i in range(pools[1])]
    ids = ["r"] + [f"w{i}" for i in range(len(word))] + ["t"]
    return {"nodes": nodes, "edges": [
        {"id": i, "source": left, "target": right, "label": None}
        for i, (left, right) in enumerate(zip(ids, ids[1:]))
    ]}


def graph_matches(pattern: dict, host: dict) -> list[dict[str, str]]:
    """Enumerate node-injective matches into one declared simple ROP host."""
    plabel = {node["id"]: node["label"] for node in pattern["nodes"]}
    hlabel = {node["id"]: node["label"] for node in host["nodes"]}
    if len(plabel) != len(pattern["nodes"]) or len(hlabel) != len(host["nodes"]):
        raise ValueError("node IDs must be unique")
    host_edges = Counter(
        (edge["source"], edge["target"], edge.get("label"))
        for edge in host["edges"]
    )
    if any(multiplicity != 1 for multiplicity in host_edges.values()):
        raise ValueError("declared ROP hosts must not contain parallel edges")
    ordered = sorted(
        plabel,
        key=lambda node: sum(label == plabel[node] for label in hlabel.values()),
    )
    found = []

    def visit(index, mapping, used):
        if index == len(ordered):
            image_edges = Counter(
                (mapping[edge["source"]], mapping[edge["target"]], edge.get("label"))
                for edge in pattern["edges"]
            )
            if not image_edges - host_edges:
                found.append(dict(mapping))
            return
        node = ordered[index]
        for target, label in hlabel.items():
            if target not in used and label == plabel[node]:
                mapping[node] = target
                used.add(target)
                visit(index + 1, mapping, used)
                used.remove(target)
                del mapping[node]

    visit(0, {}, set())
    return found


def apply_rule(rule: dict, host: dict, match: dict[str, str]) -> dict:
    """Execute one exported DPO rule interface on a declared simple ROP host."""
    nodes = {n["id"]: n["label"] for n in host["nodes"]}
    edges = Counter(
        (edge["source"], edge["target"], edge.get("label"))
        for edge in host["edges"]
    )
    if any(multiplicity != 1 for multiplicity in edges.values()):
        raise ValueError("declared ROP hosts must not contain parallel edges")
    kept_lhs = set(rule["preserved_nodes"])
    kept_edges = {item["left"] for item in rule["preserved_edges"]}
    deleted_edges = Counter(
        (match[edge["source"]], match[edge["target"]], edge.get("label"))
        for edge in rule["lhs"]["edges"]
        if edge["id"] not in kept_edges
    )
    assert not deleted_edges - edges, "matched LHS edge absent from host"
    deleted_nodes = {
        match[node]
        for node in {item["id"] for item in rule["lhs"]["nodes"]} - kept_lhs
    }
    for target in deleted_nodes:
        incident = Counter({edge: count for edge, count in edges.items() if target in edge[:2]})
        assert not incident - deleted_edges, "DPO dangling violation"
    edges.subtract(deleted_edges)
    edges += Counter()
    for target in deleted_nodes:
        del nodes[target]
    rhs_map = {right: match[left] for left, right in rule["preserved_nodes"].items()}
    serial = 0
    for node in rule["rhs"]["nodes"]:
        if node["id"] not in rhs_map:
            while f"new{serial}" in nodes:
                serial += 1
            rhs_map[node["id"]] = f"new{serial}"
            nodes[f"new{serial}"] = node["label"]
    for edge in rule["rhs"]["edges"]:
        if not any(item["right"] == edge["id"] for item in rule["preserved_edges"]):
            edges[(rhs_map[edge["source"]], rhs_map[edge["target"]], edge.get("label"))] += 1
    return {"nodes": [{"id": i, "label": label} for i, label in nodes.items()],
            "edges": [{"id": i, "source": s, "target": t, "label": label}
                      for i, (s, t, label) in enumerate(sorted(edges.elements()))]}


def decode_state(graph: dict) -> tuple[tuple[str, ...], tuple[int, int]]:
    """Validate and decode the complete one-path-plus-isolated-pools language."""
    labels = {n["id"]: n["label"] for n in graph["nodes"]}
    if len(labels) != len(graph["nodes"]):
        raise ValueError("node IDs must be unique")
    allowed_labels = {"R", "T", "A", "B", "F_A", "F_B"}
    if not set(labels.values()) <= allowed_labels:
        raise ValueError("host contains an unknown node label")
    roots = [node for node, label in labels.items() if label == "R"]
    terminals = [node for node, label in labels.items() if label == "T"]
    if len(roots) != 1 or len(terminals) != 1:
        raise ValueError("host must contain exactly one root and one terminal")

    edge_tuples = [
        (edge["source"], edge["target"], edge.get("label"))
        for edge in graph["edges"]
    ]
    if any(label is not None for _, _, label in edge_tuples):
        raise ValueError("declared ROP host edges must be unlabeled")
    if len(edge_tuples) != len(set(edge_tuples)):
        raise ValueError("declared ROP hosts must not contain parallel edges")
    if any(source not in labels or target not in labels for source, target, _ in edge_tuples):
        raise ValueError("edge endpoint absent from host")

    pool_nodes = {
        node for node, label in labels.items() if label in {"F_A", "F_B"}
    }
    if any(
        source in pool_nodes or target in pool_nodes
        for source, target, _ in edge_tuples
    ):
        raise ValueError("free-pool vertices must be isolated")
    word_nodes = set(labels) - pool_nodes
    outgoing_lists = {node: [] for node in word_nodes}
    incoming_lists = {node: [] for node in word_nodes}
    for source, target, _ in edge_tuples:
        outgoing_lists[source].append(target)
        incoming_lists[target].append(source)
    if any(len(targets) > 1 for targets in outgoing_lists.values()):
        raise ValueError("word path must not branch")
    if any(len(sources) > 1 for sources in incoming_lists.values()):
        raise ValueError("word path must not merge")

    root, terminal = roots[0], terminals[0]
    if incoming_lists[root] or outgoing_lists[terminal]:
        raise ValueError("root and terminal have invalid path degree")
    outgoing = {
        node: targets[0] for node, targets in outgoing_lists.items() if targets
    }
    current = root
    word, seen = [], set()
    while labels[current] != "T":
        if current in seen or current not in outgoing:
            raise ValueError("word component is cyclic or does not reach the terminal")
        seen.add(current)
        current = outgoing[current]
        if labels[current] != "T":
            if labels[current] not in {"A", "B"}:
                raise ValueError("internal word vertices must be A or B")
            word.append(labels[current])
    seen.add(terminal)
    if seen != word_nodes:
        raise ValueError("host contains a disconnected non-pool vertex")
    return tuple(word), (
        sum(value == "F_A" for value in labels.values()),
        sum(value == "F_B" for value in labels.values()),
    )


def occurrences(word: tuple[str, ...], motif: tuple[str, ...]) -> int:
    return sum(
        word[index : index + len(motif)] == motif
        for index in range(len(word) - len(motif) + 1)
    )


def finite_validation(inputs: dict) -> dict:
    """Execute the stated finite word/pool check, without using it as a proof."""
    assert inputs["edge_label_representation"] is None
    assert len(inputs["rules"]) == 8
    assert len({rule["name"] for rule in inputs["rules"]}) == 8
    states = [
        (tuple(word), pools)
        for length in range(4)
        for word in product("AB", repeat=length)
        for pools in product(range(3), repeat=2)
    ]
    rule_checks = disabled_checks = applications = query_checks = 0
    seeds = {seed["name"]: seed["graph"] for seed in inputs["seeds"]}
    state_results = []
    for word, pools in states:
        host = path_graph(word, pools)
        state_record = {
            "word": "".join(word),
            "pools": {"A": pools[0], "B": pools[1]},
            "rule_matches": {},
            "executed_successors": [],
            "query_embeddings": {},
        }
        for rule in inputs["rules"]:
            matches = graph_matches(rule["lhs"], host)
            symbol = "A" if rule["rate"].endswith("A") else "B"
            is_feed = rule["name"].startswith("feed_")
            previous = rule["name"].split("_")[1] if not is_feed else None
            expected_count = (
                1
                if is_feed
                else pools["AB".index(symbol)]
                * int((word[-1] if word else "R") == previous)
            )
            assert len(matches) == expected_count
            state_record["rule_matches"][rule["name"]] = len(matches)
            rule_checks += 1
            disabled_checks += int(not matches)
            decoded_targets = []
            for match in matches:
                target = decode_state(apply_rule(rule, host, match))
                expected_pools = list(pools)
                if is_feed:
                    expected_pools["AB".index(symbol)] += 1
                    expected_word = word
                else:
                    expected_pools["AB".index(symbol)] -= 1
                    expected_word = word + (symbol,)
                assert target == (expected_word, tuple(expected_pools))
                applications += 1
                decoded_targets.append(
                    {
                        "word": "".join(target[0]),
                        "pools": {"A": target[1][0], "B": target[1][1]},
                    }
                )
            if decoded_targets:
                state_record["executed_successors"].append(
                    {
                        "rule": rule["name"],
                        "multiplicity": len(decoded_targets),
                        "targets": decoded_targets,
                    }
                )
        for name in QUERIES:
            embeddings = len(graph_matches(seeds[name], host))
            assert embeddings == occurrences(word, tuple(name.split("_", 1)[1]))
            state_record["query_embeddings"][name] = embeddings
            query_checks += 1
        state_results.append(state_record)
    return {
        "status": "PASS",
        "maximum_word_length": 3,
        "pool_range_each": [0, 2],
        "states": len(states),
        "rule_match_checks": rule_checks,
        "disabled_rule_checks": disabled_checks,
        "executed_successor_applications": applications,
        "query_embedding_checks": query_checks,
        "state_results": state_results,
        "counterexamples": [],
        "scope": "finite validation only; not an all-state correspondence proof",
    }


def map_graph(graph: dict, ring):
    """Map one injectively matched graph, never a product of graph factors."""
    labels = {n["id"]: n["label"] for n in graph["nodes"]}
    free = Counter(label for label in labels.values() if label in ("F_A", "F_B"))
    nonfree = {node for node, label in labels.items() if label not in ("F_A", "F_B")}
    if any(e["source"] not in nonfree or e["target"] not in nonfree for e in graph["edges"]):
        return None, "edge incident to free-pool vertex"
    value = ring.constant(1)
    for symbol in "AB":
        for offset in range(free[f"F_{symbol}"]):
            falling_factor = ring.add(
                ring.coordinate(f"M_{symbol}"), ring.constant(-offset)
            )
            value = ring.multiply(value, falling_factor)
    if not nonfree:
        return value, None
    if any(e.get("label") is not None for e in graph["edges"]):
        return None, "unexpected edge label in unlabeled representation"
    if len({e["source"] for e in graph["edges"]}) != len(graph["edges"]):
        return None, "branching word component"
    if len({e["target"] for e in graph["edges"]}) != len(graph["edges"]):
        return None, "merging word component"
    edges = {e["source"]: e["target"] for e in graph["edges"]}
    incoming = {e["target"] for e in graph["edges"]}
    starts = [node for node in nonfree if node not in incoming]
    if len(starts) != 1 or len(graph["edges"]) != len(nonfree)-1:
        return None, "non-path or multiple word components"
    sequence = []
    node = starts[0]
    while True:
        sequence.append(labels[node])
        if node not in edges:
            break
        node = edges[node]
    if len(sequence) != len(nonfree):
        return None, "disconnected or cyclic word component"
    if set(sequence) <= {"A", "B"}:
        name = (
            f"N_{sequence[0]}"
            if len(sequence) == 1
            else "C_" + "".join(sequence)
        )
        if name not in ring.variables:
            return None, f"internal motif {name} absent from current coordinates"
        return ring.multiply(value, ring.coordinate(name)), None
    if (
        sequence[-1:] == ["T"]
        and set(sequence[:-1]) <= {"A", "B"}
        and sequence[:-1]
    ):
        name = "S_" + "".join(sequence[:-1])
        if name not in ring.variables:
            return None, f"suffix {name} absent from current coordinates"
        return ring.multiply(value, ring.coordinate(name)), None
    if sequence == ["R", "T"]:
        empty = ring.add(
            ring.constant(1),
            ring.scale(ring.coordinate("S_A"), -1),
            ring.scale(ring.coordinate("S_B"), -1),
        )
        return ring.multiply(value, empty), None
    return None, "rooted or unsupported word path"


def compare_equations(run: dict, model) -> dict:
    """Compare every equation whose graph factors all have exact mappings."""
    ring = model.ring
    channels = tuple(dict.fromkeys(event.channel for event in model.events))
    rows = []
    for row in run["equations"]:
        lhs, lhs_reason = map_graph(row["lhs"], ring)
        mapped = []
        rhs_by_channel = {channel: ring.constant(0) for channel in channels}
        reasons = []
        for term in row["terms"]:
            product_value = ring.constant(1)
            for graph in term["graphs"]:
                value, reason = map_graph(graph, ring)
                if value is None:
                    reasons.append(reason)
                    break
                product_value = ring.multiply(product_value, value)
            else:
                if len(term["rates"]) == 1:
                    channel = term["rates"][0]
                    rhs_by_channel[channel] = ring.add(
                        rhs_by_channel[channel],
                        ring.scale(product_value, str(term["coefficient"])),
                    )
                    mapped.append(True)
                    continue
                reasons.append("term does not contain exactly one rate channel")
            mapped.append(False)
        if lhs is not None and all(mapped):
            composition = {
                channel: generator_action(model, lhs, channel=channel)
                for channel in channels
            }
            equal = all(rhs_by_channel[channel] == composition[channel] for channel in channels)
            rows.append(
                {
                    "lhs": row["lhs"],
                    "mapped_lhs": ring.record(lhs),
                    "status": "EQUAL" if equal else "MISMATCH",
                    "all_rate_images_equal": equal,
                    "external_by_channel": {
                        channel: ring.record(rhs_by_channel[channel])
                        for channel in channels
                    },
                    "composition_by_channel": {
                        channel: ring.record(composition[channel])
                        for channel in channels
                    },
                }
            )
        else:
            all_reasons = ([lhs_reason] if lhs_reason else []) + reasons
            rows.append(
                {
                    "lhs": row["lhs"],
                    "status": "UNMAPPED",
                    "reasons": sorted(set(all_reasons)),
                }
            )
    equal = sum(row["status"] == "EQUAL" for row in rows)
    mismatch = sum(row["status"] == "MISMATCH" for row in rows)
    unmapped = sum(row["status"] == "UNMAPPED" for row in rows)
    status = "MISMATCH" if mismatch else ("COMPARISON_INCOMPLETE" if unmapped else "PASS")
    return {
        "status": status,
        "mapped_equal": equal,
        "mapped_mismatch": mismatch,
        "unmapped": unmapped,
        "total": len(rows),
        "rows": rows,
    }


def aggregate_comparison_status(external: list[dict], comparisons: list[dict]) -> str:
    """Aggregate every requested run; mismatch takes precedence over missing data."""
    if any(item["status"] == "MISMATCH" for item in comparisons):
        return "MISMATCH"
    if (
        len(comparisons) == len(external)
        and all(item["status"] == "PASS" for item in external)
        and all(item["status"] == "PASS" for item in comparisons)
    ):
        return "PASS"
    return "COMPARISON_INCOMPLETE"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--java", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    checkout, java, out = args.checkout.resolve(), args.java.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    fragger_identity = require_source_identity(
        checkout, expected_commit=COMMIT, snapshot_name=FRAGGER_SNAPSHOT,
        git_paths=("src/main",),
    )
    compiler = checkout / "compiler"
    jars = [compiler / f"scala-{name}-2.13.2.jar" for name in ("compiler", "library", "reflect")]
    assert {path.name: sha256(path) for path in jars} == JAR_HASHES
    classes = checkout / "sequence-comparison-classes"; classes.mkdir(exist_ok=True)
    compile_command = [java, "-Xmx1g", "-cp", str(compiler / "*"), "scala.tools.nsc.Main", "-usejavacp", "-d", classes,
                       *sorted((checkout / "src/main/scala/graph-rewriting").glob("*.scala")),
                       Path(__file__).with_name("RopSequenceHarness.scala")]
    compile_started = time.perf_counter(); compiled = command(compile_command, checkout, 120)
    compile_seconds = time.perf_counter() - compile_started
    (out / "compile.log").write_text(compiled.stdout + compiled.stderr, encoding="utf-8")
    export_command = [java, "-Xmx1g", "-cp", f"{classes}:{compiler / '*'}",
                      "hz.ricardo.graph_rewriting.RopSequenceHarness", "export"]
    exported = command(export_command, checkout, 30)
    export_lines = [line for line in exported.stdout.split("\n") if line.startswith("INPUT_JSON=")]
    if len(export_lines) != 1:
        raise ValueError("Scala harness did not emit exactly one INPUT_JSON record")
    inputs = json.loads(export_lines[0].split("=", 1)[1])
    (out / "scala-input.json").write_text(json.dumps(inputs, indent=2) + "\n", encoding="utf-8")
    external = []
    for budget in (12, 24):
        run_command = [java, "-Xmx1g", "-cp", f"{classes}:{compiler / '*'}",
                       "hz.ricardo.graph_rewriting.RopSequenceHarness", str(budget)]
        external.append(external_worklist(run_command, checkout, 60,
                                          out / f"fragger-budget-{budget}.log", budget))
    validation_started = time.perf_counter(); validation = finite_validation(inputs); validation_seconds = time.perf_counter() - validation_started
    model = projection_checks.composition_compilation("copolymer").model
    comparisons = []
    for item in external:
        if item["json"] is None:
            comparisons.append(
                {
                    "budget": item["budget"],
                    "status": "EXTERNAL_FAILURE",
                    "external_status": item["status"],
                }
            )
        else:
            comparisons.append(
                {"budget": item["budget"], **compare_equations(item["json"], model)}
            )
    comparison_status = aggregate_comparison_status(external, comparisons)
    project_identity = source_identity(ROOT)
    java_version = command([java, "-version"], checkout, 30)
    record = {"schema_version": 1, "kind": "bounded-ungated-rop-sequence-fragger-diagnostic",
              "comparison_boundary": "positive directed-word population translation; not molecular graph-input equivalence or an automatic compiler",
              "fragger": {"commit": COMMIT, "source_unmodified": True,
                           "source_identity": fragger_identity, "scala": "2.13.2",
                           "java": (java_version.stdout + java_version.stderr).strip(),
                           "jar_sha256": {path.name: sha256(path) for path in jars}},
              "project_revision": project_identity,
              "inputs": {"scala_export": inputs,
                         "language_identities": ["[R->T]=1-[A->T]-[B->T]", "same identity with disconnected free-pool-only context", "repeated F labels count by falling factorial"],
                         "transformer_supplied_to_fragger": "none; identities are documented mapping identities, not a target-basis reducer"},
              "finite_validation": validation, "external_worklists": external,
              "molecular_composition": {"input_definition_sha256": model.input_definition_sha256,
                  "all_jumps_composition_derived": model.provenance.all_jumps_composition_derived,
                  "event_classes": len(model.events), "coordinates": list(model.ring.variables)},
              "symbolic_comparison": {"status": comparison_status, "runs": comparisons,
                  "scope": "Every exactly translated graph factor is compared channel-by-channel; unsupported graph shapes remain explicit and are not zeroed or guessed."},
              "hashes": {"harness": sha256(Path(__file__).with_name("RopSequenceHarness.scala")), "runner": sha256(Path(__file__)),
                         "engine": sha256(Path(engine.__file__)), "rop_model": sha256(ROOT / "src/mfobs/models/molecular_rop/inputs.py")},
              "diagnostic_phase_seconds": {"compile": compile_seconds, "finite_validation_only": validation_seconds,
                  "external_processes": {str(item["budget"]): item["process_seconds"] for item in external}},
              "process_limits": {"compile_seconds": 120, "external_seconds": 60, "java_heap": "1g"},
              "executed_compile_command": [str(x) for x in compile_command],
              "executed_input_export_command": [str(x) for x in export_command]}
    (out / "comparison.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "DIAGNOSIS_RECORDED", "validation": validation["status"],
                      "external": [{"budget": item["budget"], "status": item["status"]} for item in external], "out": str(out)}, indent=2))


if __name__ == "__main__":
    main()
