#!/usr/bin/env python3
"""Fresh FRAGGER checks for two declared ROP population projections.

The comparison boundary is deliberately narrow.  FRAGGER alone receives the
edgeless population rules and the seed ``N``; it receives no molecular graph,
requested moment basis, or equation.  After FRAGGER has selected its left-hand
patterns, this runner independently recomputes *those supplied patterns* with
MolFoundry's edgeless DPO composition.  Separately, it derives the current
atom-resolved ROP jump model by its existing DPO-composition pipeline and
compares the mapped generator functions.  Neither MolFoundry path is a second
worklist-discovery claim here.  This is not a claim that FRAGGER natively
implements molecular NoEdge/DPO.
"""
from __future__ import annotations

import mfobs.checks.rop_projection as projection_checks

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from mfobs.core import engine
from mfobs.core.conditions import parse_observable
from mfobs.core.graph import make_graph
from mfobs.core.rules import make_rule
from mfobs.core.polynomial_jump_closure import generator_action
from mfobs.core.source_identity import (
    FRAGGER_SNAPSHOT,
    require_source_identity,
    source_identity,
)


COMMIT = "0d9595c5d6e8918fc605f7d65778568a74f5e2d4"
JAR_HASHES = {
    "scala-compiler-2.13.2.jar": "b04da83111decff5b36e87c74d73f5e2231b63ce605f2602e862bae917d7761a",
    "scala-library-2.13.2.jar": "6a1c9846feb8f5d7d33b91e61e056d272215a79747a2b7b5e5ece690bf151d5d",
    "scala-reflect-2.13.2.jar": "21a7288874f186d6049e8b14a986011ffd96cb587f206eccf946624f2d476815",
}

N, M, O, C = sp.symbols("N M O C")
N_B, M_B, X = sp.symbols("N_B M_B X")
beta_B, k_B, a, b = sp.symbols("beta_B k_B a b")
POPULATION_SPECIES = {"N": N, "M": M, "O": O, "C": C}
RATE_SYMBOLS = {str(rate): rate for rate in (beta_B, k_B, a, b)}
VARIANTS = ("baseline", "noedge")


def command(args, cwd: Path, *, timeout: int) -> tuple[str, str]:
    process = subprocess.run(
        [str(value) for value in args], cwd=cwd, text=True,
        capture_output=True, timeout=timeout, check=True,
    )
    return process.stdout, process.stderr


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def population_graph(labels: list[str] | tuple[str, ...], start: int = 0) -> engine.mf.Graph:
    return make_graph({index + start: label for index, label in enumerate(labels)})


def population_rule(name: str, left, right):
    return make_rule(name, left, right)


def population_rules(variant: str):
    common = [
        ("beta_B", population_rule("feed_B", population_graph([]), population_graph(["M"]))),
    ]
    if variant == "baseline":
        return common + [
            ("k_B", population_rule("insert_B", population_graph(["M"]), population_graph(["N"], 1))),
        ]
    if variant == "noedge":
        return common + [
            ("k_B", population_rule(
                "insert_B", population_graph(["O", "M"]), make_graph({0: "O", 2: "N"})
            )),
            ("a", population_rule("donor_bind", population_graph(["O"]), population_graph(["C"], 1))),
            ("b", population_rule("donor_unbind", population_graph(["C"]), population_graph(["O"], 1))),
        ]
    raise ValueError(f"unknown variant {variant!r}")


def row_keyed(row: dict) -> dict[tuple[tuple[str, ...], tuple[str, ...]], Fraction]:
    terms: defaultdict[tuple[tuple[str, ...], tuple[str, ...]], Fraction] = defaultdict(Fraction)
    for term in row["terms"]:
        terms[(tuple(sorted(term["rates"])), tuple(sorted(term["graph"])))] += Fraction(str(term["coefficient"]))
    return {key: value for key, value in terms.items() if value}


def composition_row(labels: list[str], rules) -> dict:
    """Obtain a raw shared-population generator image by DPO composition."""
    observable_graph = population_graph(labels)
    observable = population_rule("observable", observable_graph, observable_graph)
    coefficients: defaultdict[tuple[tuple[str, ...], tuple[str, ...]], int] = defaultdict(int)
    raw_compositions = 0
    for rate, reaction in rules:
        commutator = engine.commutator(
            observable, reaction, connected=False, include_empty=True
        )
        raw_compositions += commutator.left_raw + commutator.right_raw
        for key, coefficient in commutator.observable_coefficients.items():
            pattern, adjacency = parse_observable(commutator.observable_representatives[key])
            assert not pattern.edges and not pattern.no_edges and not adjacency
            coefficients[(rate,), tuple(sorted(pattern.labels.values()))] += coefficient
    return {
        "lhs": list(labels),
        "terms": [
            {"coefficient": value, "rates": list(rates), "graph": list(graph)}
            for (rates, graph), value in sorted(coefficients.items()) if value
        ],
        "raw_compositions": raw_compositions,
    }


def count_pattern(labels: tuple[str, ...], *, noedge: bool) -> sp.Expr:
    """Injective edgeless-pattern count, reduced only by O+C=1 when declared."""
    expression = sp.prod(
        sp.prod(POPULATION_SPECIES[label] - index for index in range(count))
        for label, count in Counter(labels).items()
    )
    expression = sp.expand(expression)
    if noedge:
        expression = expression.subs(C, 1 - O)
        # O is an indicator on the singleton-catalyst language.
        expression = expression.subs(O, 0) + O * (expression.subs(O, 1) - expression.subs(O, 0))
    return sp.expand(expression)


def row_expression(row: dict, *, noedge: bool) -> sp.Expr:
    return sp.expand(sum(
        sp.Rational(coefficient.numerator, coefficient.denominator)
        * sp.prod(RATE_SYMBOLS[rate] for rate in rates)
        * count_pattern(graph, noedge=noedge)
        for (rates, graph), coefficient in row_keyed(row).items()
    ))


def polynomial_to_sympy(value, ring, names_to_symbols: dict[str, sp.Symbol]) -> sp.Expr:
    result = sp.Integer(0)
    for monomial, coefficient in value.terms:
        term = sp.Rational(coefficient.numerator, coefficient.denominator)
        for name, power in zip(ring.variables, monomial):
            term *= names_to_symbols[name] ** power
        result += term
    return sp.expand(result)


def sympy_to_polynomial(value: sp.Expr, ring, symbols_to_names: dict[sp.Symbol, str]):
    polynomial = sp.Poly(sp.expand(value), *symbols_to_names)
    terms = {}
    for powers, coefficient in polynomial.terms():
        named_powers = {
            symbols_to_names[symbol]: exponent
            for symbol, exponent in zip(symbols_to_names, powers) if exponent
        }
        terms[tuple(
            named_powers.get(name, 0) for name in ring.variables
        )] = Fraction(str(coefficient))
    return ring.polynomial(terms)


def rop_mapping(variant: str) -> dict:
    if variant == "baseline":
        return {
            "external_to_rop": {"N": "N_B", "M": "M_B"},
            "rate_to_rop": {"beta_B": "beta_B", "k_B": "k_B"},
            "identities": [],
            "normal_forms": {
                "[N]": "N_B", "[M]": "M_B", "[empty]": "1",
            },
        }
    if variant == "noedge":
        return {
            "external_to_rop": {"N": "N_B", "M": "M_B", "O": "X", "C": "1-X"},
            "rate_to_rop": {"beta_B": "beta_B", "k_B": "k_B", "a": "a", "b": "b"},
            "identities": [
                "O+C=1", "O^2=O", "C^2=C", "O*C=0",
                "[P+C]=[P]-[P+O]", "[P+O+O]=0",
            ],
            "normal_forms": {
                "[N]": "N_B", "[M]": "M_B", "[O]": "X",
                "[C]": "1-X", "[M+O]": "M_B*X", "[empty]": "1",
            },
        }
    raise ValueError(variant)


def rop_coordinate_text(value: sp.Expr) -> str:
    """Render a common-population function in the current ROP coordinates."""
    return str(sp.expand(value.subs({N: N_B, M: M_B, O: X, C: 1 - X})))


def molecular_composition_images(variant: str, rows: list[dict]) -> dict:
    """Use the current ROP DPO-composition result, never a hand-fed jump table."""
    compiled = projection_checks.composition_compilation(variant).model
    assert compiled.provenance.all_jumps_composition_derived
    assert not compiled.provenance.supplied_delta_table
    expected_coordinates = ("N_B", "M_B") if variant == "baseline" else ("N_B", "M_B", "X")
    assert compiled.ring.variables == expected_coordinates
    expected_channels = ("beta_B", "k_B") if variant == "baseline" else ("beta_B", "k_B", "a", "b")
    assert tuple(dict.fromkeys(event.channel for event in compiled.events)) == expected_channels

    coordinate_symbols = {"N_B": N, "M_B": M, "X": O}
    symbols_to_names = {value: key for key, value in coordinate_symbols.items() if key in compiled.ring.variables}
    names_to_symbols = {key: value for value, key in symbols_to_names.items()}
    noedge = variant == "noedge"
    images = []
    for row in rows:
        observable_count = count_pattern(tuple(row["lhs"]), noedge=noedge)
        observable = sympy_to_polynomial(observable_count, compiled.ring, symbols_to_names)
        expression = sp.Integer(0)
        for channel in expected_channels:
            expression += RATE_SYMBOLS[channel] * polynomial_to_sympy(
                generator_action(compiled, observable, channel=channel),
                compiled.ring, names_to_symbols,
            )
        images.append({
            "lhs": row["lhs"],
            "mapped_observable": rop_coordinate_text(observable_count),
            "composition_generator_image_comparison_coordinates": str(sp.expand(expression)),
            "composition_generator_image_rop_coordinates": rop_coordinate_text(expression),
        })
    return {
        "input_definition_sha256": compiled.input_definition_sha256,
        "provenance": compiled.provenance.record(),
        "event_classes": [
            {
                "identifier": event.identifier, "channel": event.channel,
                "rule": event.rule, "predicate": event.predicate,
                "delta": dict(zip(compiled.ring.variables, event.delta)),
                "propensity": compiled.ring.text(event.propensity),
            }
            for event in compiled.events
        ],
        "images": images,
    }


def parse_fresh_run(checkout: Path, java: Path, classes: Path, compiler: Path, variant: str, out: Path) -> dict:
    started = time.perf_counter()
    command_line = [
        java, "-Xmx1g", "-cp", os.pathsep.join([str(classes), str(compiler / "*")]),
        "hz.ricardo.graph_rewriting.RopPopulationHarness", variant,
    ]
    stdout, stderr = command(command_line, checkout, timeout=60)
    (out / f"fragger-{variant}.log").write_text(stdout + stderr, encoding="utf-8")
    lines = [line for line in stdout.splitlines() if line.startswith("FRAGGER_JSON=")]
    assert len(lines) == 1, stdout
    run = json.loads(lines[0].split("=", 1)[1])
    run["external_process_seconds"] = time.perf_counter() - started
    run["executed_command"] = [str(value) for value in command_line]
    assert run["variant"] == variant and run["seed"] == ["N"]
    assert not run["unresolved"], run
    return run


def compare_variant(variant: str, fresh: dict) -> dict:
    """Compare actual external equations with both MolFoundry paths."""
    noedge = variant == "noedge"
    rules = population_rules(variant)
    population_rows = []
    molecular = molecular_composition_images(variant, fresh["equations"])
    comparisons = []
    for row, molecular_image in zip(fresh["equations"], molecular["images"]):
        population_row = composition_row(row["lhs"], rules)
        external_image = row_expression(row, noedge=noedge)
        population_image = row_expression(population_row, noedge=noedge)
        assert sp.expand(external_image - population_image) == 0, (row, population_row)
        molecular_image_value = sp.sympify(
            molecular_image["composition_generator_image_comparison_coordinates"],
            locals={"N": N, "M": M, "O": O, "C": C, **RATE_SYMBOLS},
        )
        assert sp.expand(external_image - molecular_image_value) == 0, (row, molecular_image)
        population_rows.append(population_row)
        comparisons.append({
            "lhs": row["lhs"],
            "mapped_observable": molecular_image["mapped_observable"],
            "fresh_fragger_equals_population_dpo_composition": True,
            "fresh_fragger_equals_molecular_rop_composition_generator": True,
            "population_dpo_raw_compositions": population_row["raw_compositions"],
            "fresh_equation": str(external_image),
            "molecular_generator_image": molecular_image["composition_generator_image_rop_coordinates"],
        })
    expected_dimension = 3 if variant == "baseline" else 5
    assert fresh["ode_count"] == expected_dimension
    return {
        "variant": variant,
        "population_projection": {
            "state_language": (
                "finite edgeless N/M populations" if not noedge else
                "finite edgeless N/M/O/C populations satisfying O+C=1"
            ),
            "rules": [
                {"rate": rate, "left_right_rule": reaction.getGMLString()}
                for rate, reaction in rules
            ],
            "seed": ["N"],
            "per_match_normalization": "injective labelled embeddings; no automorphism division; one per-match multiplier",
            "mapping_to_current_rop": rop_mapping(variant),
        },
        "fresh_fragger": fresh,
        "shared_population_molfoundry_dpo_composition": population_rows,
        "current_molecular_rop_composition": molecular,
        "comparable_generator_images": comparisons,
        "closure": {
            "status": "CLOSED",
            "unresolved_external_observables": fresh["unresolved"],
            "basis_dimension_including_constant": expected_dimension,
            "basis_mapping": [
                {"external_pattern": row["lhs"], "rop_function": image["mapped_observable"]}
                for row, image in zip(fresh["equations"], molecular["images"])
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--java", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    checkout, java, out = args.checkout.resolve(), args.java.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    compile_started = time.perf_counter()
    fragger_identity = require_source_identity(
        checkout, expected_commit=COMMIT, snapshot_name=FRAGGER_SNAPSHOT,
        git_paths=("src/main",),
    )
    commit = fragger_identity["commit"]
    compiler = checkout / "compiler"
    jars = [compiler / f"scala-{part}-2.13.2.jar" for part in ("compiler", "library", "reflect")]
    assert all(path.is_file() for path in jars), "provide pinned Scala 2.13.2 jars under CHECKOUT/compiler"
    assert {path.name: sha256(path) for path in jars} == JAR_HASHES, "Scala jar checksum mismatch"
    classes = checkout / "rop-population-classes"
    classes.mkdir(exist_ok=True)
    compile_command = [
        java, "-Xmx1g", "-cp", str(compiler / "*"), "scala.tools.nsc.Main", "-usejavacp",
        "-d", classes,
        *sorted((checkout / "src/main/scala/graph-rewriting").glob("*.scala")),
        Path(__file__).with_name("RopPopulationHarness.scala"),
    ]
    compiled = command(compile_command, checkout, timeout=120)
    (out / "compile.log").write_text(compiled[0] + compiled[1], encoding="utf-8")
    compilation_seconds = time.perf_counter() - compile_started

    fresh_runs = [parse_fresh_run(checkout, java, classes, compiler, variant, out) for variant in VARIANTS]
    comparison_started = time.perf_counter()
    comparisons = [compare_variant(variant, fresh) for variant, fresh in zip(VARIANTS, fresh_runs)]
    comparison_seconds = time.perf_counter() - comparison_started
    project_identity = source_identity(ROOT)
    metadata = {
        "schema_version": 1,
        "kind": "fresh-fragger-rop-population-projection-comparison",
        "comparison_boundary": (
            "Shared inputs are reasoned edgeless population projections. Equality is asserted only for their "
            "mapped generator images; it is not molecular graph-input equivalence and does not transfer native NoEdge/NAC semantics."
        ),
        "fragger": {
            "repository": "https://github.com/rhz/graph-rewriting",
            "commit": commit,
            "source_unmodified": True,
            "source_identity": fragger_identity,
            "scala_version": "2.13.2",
            "compiler_jar_sha256": {path.name: sha256(path) for path in jars},
            "java_version": "".join(command([java, "-version"], checkout, timeout=30)).strip(),
        },
        "molfoundry": {
            "version": importlib.metadata.version("molfoundry"),
            "engine_source_sha256": sha256(Path(engine.__file__)),
            "javascript_sha256": sha256(Path(engine.mf.__file__).parent / "bin/molfoundry.js"),
            "rop_model_source_sha256": sha256(ROOT / "src/mfobs/models/molecular_rop/inputs.py"),
        },
        "project_revision": project_identity,
        "harness_sha256": sha256(Path(__file__).with_name("RopPopulationHarness.scala")),
        "runner_sha256": sha256(Path(__file__)),
        "environment": {"python": sys.version, "platform": platform.platform()},
        "diagnostic_phase_seconds": {
            "compilation": compilation_seconds,
            "fresh_external_process_wall_clock": {
                run["variant"]: run["external_process_seconds"] for run in fresh_runs
            },
            "fresh_fragger_internal_algorithm": {
                run["variant"]: run["seconds"] for run in fresh_runs
            },
            "comparison_analysis_python_and_molfoundry_only": comparison_seconds,
        },
        "timing_scope": "diagnostic phases only; no cross-tool performance or speedup claim",
        "subprocess_limits": {"compile_seconds": 120, "fresh_run_seconds": 60, "java_heap": "1g"},
        "executed_commands": {
            "compile": [str(value) for value in compile_command],
            "fresh_fragger": {run["variant"]: run["executed_command"] for run in fresh_runs},
        },
        "variants": comparisons,
    }
    (out / "comparison.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS", "variants": list(VARIANTS),
        "fragger_commit": commit,
        "basis_dimensions_including_constant": {
            row["variant"]: row["closure"]["basis_dimension_including_constant"] for row in comparisons
        },
        "out": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
