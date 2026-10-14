"""Public application adapters stay thin and route through the shared core models."""

from mfobs import (
    BrusselatorROP,
    EightSlots,
    ABA,
    RemyTrees,
    analyze,
)


def test_remy_public_adapter_delegates(monkeypatch):
    from mfobs.models.remy import analysis as remy_trees

    seen = {}

    def fake_common_analysis(**kwargs):
        seen.update(kwargs)
        return {"linear": {"status": "CLOSED"}}

    monkeypatch.setattr(remy_trees, "common_analysis", fake_common_analysis)

    backend = RemyTrees(4)
    result = analyze(backend.specification(), backend=backend, linear_limit=9, record_patterns=True).record()

    assert seen == {
        "linear_limit": 9,
        "order": 4,
        "record_patterns": True,
        "state_space": backend.specification().omega,
    }
    assert result["observable_designation"]["selected"] == "P4"


def test_aba_public_adapter_delegates(monkeypatch):
    from mfobs.models.aba import analysis as aba_analysis

    seen = {}

    def fake_analyze(selection, **kwargs):
        seen["selection"] = selection
        seen.update(kwargs)
        return {"pjc": {"status": "CLOSED"}}

    monkeypatch.setattr(aba_analysis, "analyze", fake_analyze)

    backend = ABA(('B',), alphabet=('B',), gated=False)
    result = analyze(backend.specification(), backend=backend, linear_limit=7, diagnostic_timings=True).record()

    assert seen["selection"].signature() == backend.specification().signature()
    assert seen["linear_limit"] == 7
    assert seen["diagnostic_timings"] is True
    assert result["pjc"]["status"] == "CLOSED"


def test_brusselator_public_adapter_delegates(monkeypatch):
    from mfobs.models import brusselator_pjc

    seen = {}

    def fake_analyze(**kwargs):
        seen.update(kwargs)
        return {"pjc": {"status": "CLOSED"}}

    monkeypatch.setattr(brusselator_pjc, "analyze", fake_analyze)

    backend = BrusselatorROP(("X", "Y"))
    result = analyze(backend.specification(), backend=backend, linear_limit=11).record()

    assert seen == {
        "linear_limit": 11,
        "selected_queries": ("X", "Y"),
        "state_space": backend.specification().omega,
    }
    assert result["pjc"]["status"] == "CLOSED"


def test_eight_slots_public_adapter_uses_shared_worklist(monkeypatch):
    from mfobs.core import worklist
    from mfobs import inputs

    seen = {}

    from pathlib import Path
    model = inputs.load(Path(__file__).parents[1]/"examples/eight-slots/model.json")
    monkeypatch.setattr(inputs, "load", lambda path: model)

    def fake_run(model, **kwargs):
        seen["model"] = model
        seen.update(kwargs)
        return {"status": "CLOSED"}

    monkeypatch.setattr(worklist, "run", fake_run)

    backend = EightSlots()
    result = analyze(backend.specification(), backend=backend, max_expansions=12).record()

    assert result["status"] == "CLOSED"
    assert result["source_sha256"]
    assert seen["model"] is model
    assert seen["max_expansions"] == 12
    assert seen["progress"] is False


def test_rop_public_application_exposes_the_same_inputs_used_for_analysis():
    backend = ABA.motif('ABA', gated=True)

    assert backend.motifs == ("ABA",)
    assert backend.gated is True
    assert "insert_A" in backend.rule_gmls()
    assert "donor_bind" in backend.rule_gmls()
    assert set(backend.query_gmls()) == {"ABA"}
    assert set(ABA(('B',), alphabet=('B',), gated=False).query_gmls()) == {"B"}


def test_brusselator_public_application_exposes_artifact_inputs():
    backend = BrusselatorROP()

    assert backend.contract().rules
    assert set(backend.query_gmls()) == {"X", "Y"}
    assert set(backend.query_gmls()).issubset(backend.coordinate_gmls())


def test_path_and_molecular_models_have_explicitly_different_domains():
    from mfobs import MolecularROP
    path = ABA.motif("XYX", alphabet=("X", "Y"))
    assert path.specification().omega.kind == "terminal-paths"
    molecule = MolecularROP.motif("ABA")
    assert molecule.specification().omega.kind == "molecular-paths"
    import pytest
    with pytest.raises(ValueError, match="chemical repeat"):
        MolecularROP.motif("XYX", alphabet=("X", "Y"))
