"""A new graph family must work without model-specific synthesis Python."""
import importlib.util
from pathlib import Path
import shutil

import pytest
from mfobs.inputs import load
from mfobs.core.worklist import run

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "examples/gated-path"
spec = importlib.util.spec_from_file_location("gated_path_independent_checks", ROOT / "scripts/validation/verify_gated_path.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@pytest.fixture(scope="module")
def gated_result():
    return run(load(SOURCE / "model.json"))


def test_new_declarative_family_discovers_from_one_seed(gated_result):
    result = gated_result
    assert result["closed"] and not result["unknown"]
    assert len(result["basis"]) == 3
    assert result["basis_counting_functions"]["V"] == {"power|0|0|0": "1", "power|1|0|0": "1"}
    report = checker.verify(load(SOURCE / "model.json"), result, max_n=5)
    assert report["symbolic_generator_identities_all_n"] == 9
    assert report["independent_graph_generator_identities"] == 216


def test_noedge_is_dynamically_active_not_redundant(tmp_path, gated_result):
    destination = tmp_path / "ungated"
    shutil.copytree(SOURCE, destination)
    rule = destination / "rules/subdivide.gml"
    rule.write_text(rule.read_text(encoding="utf-8").replace("  constrainNoEdge [ source 3 target 5 ]\n", ""), encoding="utf-8")
    inputs = load(destination / "model.json")
    result = run(inputs)
    assert result["closed"] and len(result["basis"]) == 2
    assert len(result["basis"]) != len(gated_result["basis"])
    assert all(key.endswith("|0") for form in result["basis_counting_functions"].values() for key in form)
    report = checker.verify(inputs, result, max_n=4, gated=False)
    assert report["passed"]


def test_compiled_constraints_match_python_predicate_counts():
    from mfobs.core import engine
    from mfobs.core.composition import group_compositions
    from mfobs.core.forbidden import LocalPatternPredicate
    inputs = load(SOURCE / 'model.json')
    observable = inputs.engine_observables()['V']
    specification = inputs.forbidden_spec().engine_payload()
    predicate = LocalPatternPredicate(specification)
    for reaction in inputs.engine_rules().values():
        for first, second in ((observable, reaction), (reaction, observable)):
            fast = engine.grouped_restricted_common_compositions(
                first, second, forbidden_specification=specification)
            public = group_compositions(first, second, violation=predicate.violation)
            assert fast[1] == public[1]
            assert fast[0].accepted_raw == public[0].accepted_raw


@pytest.mark.parametrize("degree", [2, 4])
def test_higher_moments_discovered_without_basis_input(degree):
    inputs = load(SOURCE / f"model-degree{degree}.json")
    result = run(inputs)
    assert result["closed"] and len(result["basis"]) == checker.analytic_cyclic_rank(degree)["rank"]
    assert len(inputs.data["observables"]) == len(inputs.data["initial_observables"]) == 1
    assert checker.verify(inputs, result, max_n=4)["passed"]


def test_nonclosing_feedback_returns_unknown_not_truncated_ode():
    inputs = load(SOURCE / "model-feedback.json")
    result = run(inputs, max_expansions=5)
    assert not result["closed"] and result["status"] == "UNKNOWN"
    assert result["matrix"] is None and len(result["basis"]) == 8
    assert result["pending"]
    assert checker.verify(inputs, result, max_n=4, feedback=True, partial=True)["passed"]


def test_feedback_hierarchy_leading_terms_independent_of_synthesis():
    # Regression of the general induction given in docs/reusable-family.md.
    # This finite check supplements, and does not substitute for, that proof.
    sp, n, x = checker.sp, checker.n, checker.x
    lam, a, b = sp.symbols("lambda_g a b")
    f = n + 1
    for k in range(1, 9):
        f = sp.expand(lam * checker.oracle(f, "subdivide", True, True)
                      + a * checker.oracle(f, "close", True, True)
                      + b * checker.oracle(f, "open", True, True))
        polynomial = sp.Poly(f, n)
        assert polynomial.degree() == k
        assert polynomial.coeff_monomial(n**k) == lam * (-a)**(k - 1) * x


def test_generated_graph_basis_is_not_promised_to_be_minimal():
    for degree in (2, 4, 10, 11):
        diagnostic = checker.analytic_cyclic_rank(degree)
        assert diagnostic["rank"] == 2 * (degree - 1)
        assert diagnostic["rank"] < diagnostic["polynomial_space_dimension"]
