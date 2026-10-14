"""Native declarations and family execution have independent lifecycles."""
from dataclasses import fields
import pytest
import molfoundry as mf

from mfobs import (AnalysisResult, Chlorination, ModelInputs, Observable, RatedRule,
                   SimpleGraphs, SimpleGraphBackend, analyze)
from mfobs.inputs import load, ModelManifest
from mfobs.backends import ManifestBackend


def test_inputs_can_be_constructed_without_an_application_or_execution_method():
    omega = SimpleGraphs(vertex_labels={"A"}, edge_labels={"-"})
    rule = mf.ruleGMLString('rule [ruleID "birth" right[node[id 0 label "A"]]]', add=False)
    query = mf.graphGMLString('graph[node[id 4 label "A"]]', add=False)
    inputs = ModelInputs(omega, {"birth": RatedRule(rule)}, {"A": Observable("A", query)})
    assert {field.name for field in fields(inputs)} == {"omega", "rules", "observables"}
    assert not hasattr(omega, "analyze") and not hasattr(inputs, "analyze")
    result = analyze(inputs, backend=SimpleGraphBackend(omega), uniform_matches=True)
    assert isinstance(result, AnalysisResult)
    assert result.inputs is inputs
    assert result.linear.status == "CLOSED"
    assert result.support.status == "POLYNOMIAL_SUPPORT_CLOSED"
    assert result.pjc.status == "CLOSED"


def test_execution_requires_an_explicit_backend():
    backend = Chlorination()
    with pytest.raises(TypeError, match="backend"):
        analyze(backend.specification())
    with pytest.raises(TypeError, match="ModelInputs"):
        analyze(backend, backend=backend)


def test_backend_receives_the_actual_supplied_native_objects():
    family = Chlorination()
    inputs = family.specification()
    class Recorder:
        def run(self, supplied):
            assert supplied is inputs
            assert supplied.rules["chlorinate"].rule is inputs.rules["chlorinate"].rule
            assert supplied.observables["Z"] is inputs.observables["Z"]
            return {"linear": {"status": "UNKNOWN"}, "pjc": {"status": "NOT_RUN"}}
    result = analyze(inputs, backend=Recorder())
    assert result.linear.status == "UNKNOWN"
    assert result.support is None
    assert result.pjc.status == "NOT_RUN"


def test_manifest_loads_data_without_becoming_an_analyzer():
    from pathlib import Path
    manifest = load(Path(__file__).parents[1]/"examples/tethered-nac/model-guarded.json")
    assert isinstance(manifest, ModelManifest)
    assert not hasattr(manifest, "analyze")
    result = analyze(manifest.specification(), backend=ManifestBackend(manifest), max_expansions=4)
    assert result.linear.status == "CLOSED_MODULO_FINITE_COMPONENT_MULTISET_LANGUAGE"
    assert result.support is None and result.pjc is None
    assert result.linear.data["matrix"] == [["0", "-2*k"], ["0", "-k"]]
