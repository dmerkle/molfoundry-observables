"""Mutation-sensitive checks for the independent M3b graph executor."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from scripts.comparison import run_rop_sequence_comparison as comparison
from scripts.comparison.run_rop_sequence_comparison import apply_rule, decode_state, graph_matches, path_graph


def append_rule(previous="B", symbol="A"):
    return {
        "rate": f"k_{symbol}",
        "lhs": {
            "nodes": [
                {"id": "c", "label": previous},
                {"id": "t", "label": "T"},
                {"id": "f", "label": f"F_{symbol}"},
            ],
            "edges": [{"id": 1, "source": "c", "target": "t", "label": None}],
        },
        "rhs": {
            "nodes": [
                {"id": "c", "label": previous},
                {"id": "x", "label": symbol},
                {"id": "t", "label": "T"},
            ],
            "edges": [
                {"id": 2, "source": "c", "target": "x", "label": None},
                {"id": 3, "source": "x", "target": "t", "label": None},
            ],
        },
        "preserved_nodes": {"c": "c", "t": "t"},
        "preserved_edges": [],
    }


def feed_rule(symbol="A"):
    return {
        "rate": f"beta_{symbol}",
        "lhs": {"nodes": [], "edges": []},
        "rhs": {
            "nodes": [{"id": "f", "label": f"F_{symbol}"}],
            "edges": [],
        },
        "preserved_nodes": {},
        "preserved_edges": [],
    }


def test_pool_two_yields_two_matches_and_real_successors():
    host, rule = path_graph(("A", "B"), (2, 0)), append_rule()
    matches = graph_matches(rule["lhs"], host)
    assert len(matches) == 2
    assert [decode_state(apply_rule(rule, host, match)) for match in matches] == [
        (("A", "B", "A"), (1, 0))
    ] * 2


def test_wrong_orientation_symbol_and_interface_fail():
    host, rule = path_graph(("A", "B"), (1, 0)), append_rule()
    reversed_lhs = {
        **rule["lhs"],
        "edges": [{"id": 1, "source": "t", "target": "c", "label": None}],
    }
    assert graph_matches(reversed_lhs, host) == []
    assert graph_matches(append_rule(symbol="B")["lhs"], host) == []
    broken = {**rule, "preserved_nodes": {"c": "t", "t": "c"}}
    try:
        decode_state(apply_rule(broken, host, graph_matches(rule["lhs"], host)[0]))
    except (AssertionError, KeyError, ValueError):
        pass
    else:
        raise AssertionError("wrong interface was not detected")


def test_query_embedding_is_oriented():
    host = path_graph(("A", "B", "A"), (0, 0))
    query = {
        "nodes": [
            {"id": "a0", "label": "A"},
            {"id": "b", "label": "B"},
            {"id": "a1", "label": "A"},
        ],
        "edges": [
            {"id": 0, "source": "a0", "target": "b", "label": None},
            {"id": 1, "source": "b", "target": "a1", "label": None},
        ],
    }
    assert len(graph_matches(query, host)) == 1
    query["edges"][0] = {
        "id": 0,
        "source": "b",
        "target": "a0",
        "label": None,
    }
    assert graph_matches(query, host) == []


def test_matcher_rejects_parallel_host_edges_and_cannot_match_parallel_pattern():
    nodes = [{"id": "a", "label": "A"}, {"id": "b", "label": "B"}]
    one_edge = {"nodes": nodes, "edges": [{"source": "a", "target": "b", "label": None}]}
    two_edges = {"nodes": nodes, "edges": one_edge["edges"] * 2}
    assert len(graph_matches(one_edge, one_edge)) == 1
    assert graph_matches(two_edges, one_edge) == []
    with pytest.raises(ValueError, match="parallel"):
        graph_matches(one_edge, two_edges)


def test_decoder_rejects_mutated_feed_rhs_and_extraneous_host_structure():
    host = path_graph((), (0, 0))
    rule = feed_rule()
    match = graph_matches(rule["lhs"], host)[0]
    assert decode_state(apply_rule(rule, host, match)) == ((), (1, 0))

    mutated_rule = {
        **rule,
        "rhs": {
            **rule["rhs"],
            "nodes": rule["rhs"]["nodes"]
            + [{"id": "extra_orphan", "label": "A"}],
        },
    }
    mutated = apply_rule(mutated_rule, host, match)
    with pytest.raises(ValueError, match="disconnected"):
        decode_state(mutated)

    extra_anchor = path_graph(("A",), (0, 0))
    extra_anchor["nodes"].append({"id": "extra_root", "label": "R"})
    with pytest.raises(ValueError, match="exactly one root"):
        decode_state(extra_anchor)

    pool_edge = path_graph(("A",), (1, 0))
    pool_edge["edges"].append(
        {"id": 99, "source": "fa0", "target": "t", "label": None}
    )
    with pytest.raises(ValueError, match="isolated"):
        decode_state(pool_edge)


def test_external_result_distinguishes_unicode_json_parse_and_timeout(monkeypatch, tmp_path):
    payload = {"equations": [], "node": "contains\u2028separator"}
    monkeypatch.setattr(
        comparison,
        "command",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            [], 0, f"FRAGGER_JSON={json.dumps(payload, ensure_ascii=False)}\n", ""
        ),
    )
    parsed = comparison.external_worklist(
        ["fake"], tmp_path, 1, tmp_path / "parsed.log", 12
    )
    assert parsed["status"] == "PASS" and parsed["json"] == payload

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("fake", 1, output="partial", stderr="diagnostic")

    monkeypatch.setattr(comparison, "command", timeout)
    timed = comparison.external_worklist(
        ["fake"], tmp_path, 1, tmp_path / "timeout.log", 24
    )
    assert timed["status"] == "TIMEOUT" and timed["returncode"] is None
    assert (tmp_path / "timeout.log").read_text(encoding="utf-8") == "partialdiagnostic"


def test_comparison_aggregation_requires_every_requested_run():
    for failure in ("TIMEOUT", "PARSE_ERROR"):
        external = [
            {"budget": 12, "status": "PASS"},
            {"budget": 24, "status": failure},
        ]
        comparisons = [
            {"budget": 12, "status": "PASS"},
            {"budget": 24, "status": "EXTERNAL_FAILURE"},
        ]
        assert (
            comparison.aggregate_comparison_status(external, comparisons)
            == "COMPARISON_INCOMPLETE"
        )

    external[1]["status"] = "PASS"
    comparisons[1]["status"] = "PASS"
    assert comparison.aggregate_comparison_status(external, comparisons) == "PASS"

    comparisons[1]["status"] = "MISMATCH"
    assert comparison.aggregate_comparison_status(external, comparisons) == "MISMATCH"
