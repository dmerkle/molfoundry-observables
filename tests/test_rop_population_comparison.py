"""Fast checks for the ROP shared-population comparison adapter.

These tests intentionally do not replace the fresh upstream execution in the
runner.  They guard the declared projection, normalization, and functional
mapping used to compare the separate graph languages.
"""
import sympy as sp

from scripts.comparison.run_rop_population_comparison import (
    M,
    O,
    composition_row,
    count_pattern,
    molecular_composition_images,
    population_rules,
    rop_mapping,
    row_expression,
)


def test_noedge_singleton_normal_forms_are_functional_not_factorized():
    assert count_pattern(("O", "O", "M"), noedge=True) == 0
    assert count_pattern(("C", "M"), noedge=True) == M - M * O
    assert count_pattern(("O", "M"), noedge=True) == M * O


def test_baseline_population_composition_has_only_current_channels():
    row = composition_row(["N"], population_rules("baseline"))
    assert str(row_expression(row, noedge=False)) == "M*k_B"
    assert rop_mapping("baseline")["external_to_rop"] == {"N": "N_B", "M": "M_B"}


def test_noedge_population_composition_retains_the_joint_open_pool_count():
    row = composition_row(["N"], population_rules("noedge"))
    assert row_expression(row, noedge=True) == M * O * sp.Symbol("k_B")
    assert rop_mapping("noedge")["normal_forms"]["[M+O]"] == "M_B*X"


def test_molecular_pipeline_is_composition_derived_and_uses_exact_channels():
    rows = [{"lhs": ["N"], "terms": []}]
    for variant, expected in (("baseline", ("beta_B", "k_B")), ("noedge", ("beta_B", "k_B", "a", "b"))):
        result = molecular_composition_images(variant, rows)
        assert result["provenance"]["all_jumps_composition_derived"]
        assert tuple(dict.fromkeys(event["channel"] for event in result["event_classes"])) == expected
