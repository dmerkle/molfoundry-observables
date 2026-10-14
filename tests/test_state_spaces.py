"""Canonical Ω values feed model constructors and restriction checks."""

from mfobs.core.graph import make_graph
import mfobs.models.molecular_rop.inputs as molecular_inputs
import mfobs.models.remy.analysis as remy_analysis
import mfobs.models.dna_walker.analysis as walker_analysis
from dataclasses import replace
from pathlib import Path

import molfoundry as mf
import pytest

from mfobs.models.brusselator_growth_chemistry import state_family as brusselator
from mfobs.models.chlorination.inputs import state_family as chlorinated_rings
from mfobs.state_space import covalent_compositions
from mfobs.models.dna_walker.inputs import state_family as dna_walker
from mfobs.state_space import from_manifest
from mfobs.models.molecular_rop.inputs import state_family as molecular_paths
from mfobs.models.remy.inputs import state_family as planar_binary_trees
from mfobs.models.aba.inputs import state_family as terminal_paths

from mfobs.inputs import load

ROOT = Path(__file__).resolve().parents[1]


def changed(omega, **fields):
    return replace(omega, definition={**omega.record()["definition"], **fields})


def test_declarations_are_immutable_and_records_detached():
    omega = molecular_paths()
    with pytest.raises(TypeError):
        omega.definition["repeats"]["A"]["spacer_carbons"] = 7
    record = omega.record()
    record["definition"]["repeats"]["A"]["spacer_carbons"] = 7
    record["chemistry"]["constraints"]["maxValence"]["C"] = 100
    assert omega.definition["repeats"]["A"]["spacer_carbons"] == 4
    assert omega.chemical_constraints()["maxValence"]["C"] == 8


def test_shared_chemistry_validates_native_graphs_and_both_rule_sides():
    omega = covalent_compositions()
    omega.validate_graph(mf.Graph.fromGMLString(
        'graph [ node [ id 0 label "C" ] node [ id 1 label "O" ] '
        'edge [ source 0 target 1 label "=" ] ]', add=False))
    bad = mf.Rule.fromGMLString(
        'rule [ context [ node [ id 0 label "H" ] node [ id 1 label "H" ] ] '
        'right [ edge [ source 0 target 1 label "=" ] ] ]', add=False)
    with pytest.raises(ValueError, match="forbidden-bond"):
        omega.validate_rule(bad)
    over = make_graph({0: "C", **{i: "H" for i in range(1, 6)}},
                      [(0, i, "-") for i in range(1, 6)])
    with pytest.raises(ValueError, match="forbidden-valence"):
        omega.validate_graph(over)


def test_charged_aromatic_constructor_uses_shared_upper_bounds():
    from mfobs.models.chlorination.inputs import molecule
    omega = chlorinated_rings()
    for n, mask in [(0, 0), (2, 15), (4, 5)]:
        omega.validate_graph(molecule(n, mask))
    assert omega.chemistry["source"] == brusselator().chemistry["source"]


def test_coordination_and_boundary_markers_are_explicit_extensions():
    from mfobs.models.molecular_rop import inputs as rop
    omega = molecular_paths(gated=True)
    omega.validate_graph(rop.catalyst(True))
    omega.validate_graph(rop.reservoir("A"))
    with pytest.raises(ValueError, match="unknown-atom-label"):
        covalent_compositions().validate_graph(rop.catalyst(True))
    with pytest.raises(ValueError, match="forbidden-bond"):
        omega.validate_graph(make_graph({0: "Feed-A", 1: "C"}, [(0, 1, "-")]))


def test_path_alphabet_and_gate_drive_literal_rules():
    from mfobs.models.aba.inputs import source_rules
    omega = terminal_paths(("A", "C"), True)
    rules = source_rules(("A", "C"), True, state_space=omega)
    assert [r.identifier for r in rules] == [
        "feed_A", "feed_C", "insert_A", "insert_C", "donor_bind", "donor_unbind"]
    assert 'Symbol:C' in rules[3].rule.getGMLString()
    with pytest.raises(ValueError, match="unsupported state space"):
        source_rules(state_space=changed(terminal_paths(), end="two ends"))


def test_molecular_constructor_and_pruning_consume_the_declared_family():
    from mfobs.models.molecular_rop import inputs as rop
    from mfobs.models.molecular_rop import counting as molecular
    selected = rop.configured_variant(motifs=("ABA",), gated=True)
    omega = molecular_paths(selected.alphabet, selected.gated)
    for source in molecular_inputs.source_rules(selected, state_space=omega):
        omega.validate_rule(source.rule)
    assert rop.SPACER_CARBONS == {
        a: v["spacer_carbons"] for a, v in omega.definition["repeats"].items()}
    constraints = molecular_inputs.forbidden_specification(selected)
    assert constraints["maxValence"] == omega.chemical_constraints()["maxValence"]
    assert "allowedIncidentProfiles" in constraints
    with pytest.raises(ValueError, match="unsupported state space"):
        molecular.MolecularReducer(selected, state_space=changed(omega, path="two chains"))


def test_chlorine_inventory_is_used_and_wrong_topology_rejected():
    from mfobs.models.chlorination.analysis import Chlorination
    from mfobs.models.chlorination.counting import LocalConditionalReducer
    omega = chlorinated_rings(5)
    assert LocalConditionalReducer(5, state_space=omega).chlorine_total == 5
    with pytest.raises(ValueError, match="fixes total chlorine"):
        Chlorination(chlorine_total=5, state_space=omega)
    with pytest.raises(ValueError, match="unsupported state space"):
        LocalConditionalReducer(state_space=changed(chlorinated_rings(), ring_vertices=6))


def test_remy_and_dna_use_their_canonical_restrictions():
    from mfobs.models.remy import inputs as remy_trees
    from mfobs.models.dna_walker import counting as walker
    assert remy_trees.TYPED_INCIDENCE_SPECIFICATION == planar_binary_trees().record()["definition"]["incidence"]
    assert walker.MINIMUM_CYCLE_LENGTH == dna_walker().definition["minimum_cycle_length"]
    with pytest.raises(ValueError, match="unsupported state space"):
        walker_analysis.analyze([], {}, state_space=changed(dna_walker(), minimum_cycle_length=4))
    with pytest.raises(ValueError, match="unsupported state space"):
        remy_analysis.common_analysis(state_space=changed(planar_binary_trees(), children=["L"]))


def test_brusselator_projection_consumes_population_declaration():
    from mfobs.models.brusselator_growth_projection import BrusselatorGrowthGraphModel
    from mfobs.models.brusselator_growth_chemistry import rules
    omega = brusselator()
    model = BrusselatorGrowthGraphModel(input_definition_sha256="test", state_space=omega)
    assert model.omega is omega
    assert set(model.molecular_catalogue()) == set(omega.definition["molecules"])
    for rule in rules().values():
        omega.validate_rule(rule)
    with pytest.raises(ValueError, match="unsupported state space"):
        BrusselatorGrowthGraphModel(input_definition_sha256="test",
            state_space=changed(omega, persistent_components=[]))


@pytest.mark.parametrize("name", ["eight-slots", "ring8", "structural-nac-oscillator"])
def test_manifest_file_remains_the_state_language(name):
    inputs = load(ROOT / "examples" / name / "model.json")
    omega = from_manifest(inputs)
    assert omega.record()["definition"] == inputs.read("state_language")
    for rule in inputs.engine_rules().values():
        omega.validate_rule(rule)


def test_manifest_chemistry_reference_is_not_ignored(tmp_path):
    import json
    from mfobs.inputs import ModelManifest
    source = ROOT / "src/mfobs/data/ordinary-covalent-v1.json"
    chemistry = json.loads(source.read_text(encoding="utf-8"))
    chemistry["atomTypes"]["C"]["maxValence"] = "3"
    (tmp_path / "chemistry.json").write_text(json.dumps(chemistry), encoding="utf-8")
    (tmp_path / "domain.json").write_text('{"kind":"custom-local-family"}', encoding="utf-8")
    inputs = ModelManifest(tmp_path / "model.json", {
        "id": "custom", "adapter": "custom:Adapter", "state_language": "domain.json",
        "specifications": {"chemistry": "chemistry.json", "domain": "domain.json"}})
    omega = from_manifest(inputs)
    methane = make_graph({0: "C", **{i: "H" for i in range(1, 5)}},
                         [(0, i, "-") for i in range(1, 5)])
    with pytest.raises(ValueError, match="forbidden-valence"):
        omega.validate_graph(methane)


def test_all_maintained_manifests_expose_native_inputs_and_validate_chemistry():
    import json
    count = 0
    for path in sorted((ROOT / "examples").rglob("*model*.json")):
        if "adapter" not in json.loads(path.read_text(encoding="utf-8")):
            continue
        manifest = load(path)
        inputs = manifest.specification()
        assert inputs.omega.record()["definition"] == manifest.read("state_language")
        assert list(inputs.observables) == manifest.data["initial_observables"]
        for declaration in inputs.rules.values():
            assert isinstance(declaration.rule, mf.Rule)
            inputs.omega.validate_rule(declaration.rule)
        for query in inputs.observables.values():
            for coefficient, rule in query.identity_terms():
                assert coefficient == 1  # Manifest normalization remains output/readout metadata.
                inputs.omega.validate_rule(rule)
        count += 1
    assert count, "No manifest inputs were checked"


def test_manifest_specification_runs_exact_tethered_worklist():
    from mfobs import analyze, ManifestBackend
    manifest = load(ROOT / "examples/tethered-nac/model-guarded.json")
    inputs = manifest.specification()
    result = analyze(inputs, backend=ManifestBackend(manifest), max_expansions=4).record()
    assert result["closed"]
    assert result["matrix"] == [["0", "-2*k"], ["0", "-k"]]
    assert result["initial_observables"] == ["B"]
    with pytest.raises(ValueError, match="unsupported"):
        analyze(replace(inputs, omega=changed(inputs.omega, components={})), backend=ManifestBackend(manifest), max_expansions=0)


def test_manifest_specification_preserves_explicit_rate_factor():
    from fractions import Fraction
    from mfobs import analyze, ManifestBackend
    manifest = load(ROOT / "examples/tethered-nac/model-guarded.json")
    inputs = manifest.specification()
    assert inputs.rules["cyclize"].per_match == Fraction(1, 2)
    rules = dict(inputs.rules)
    rules["cyclize"] = replace(rules["cyclize"], per_match=1)
    with pytest.raises(ValueError, match="unsupported"):
        analyze(replace(inputs, rules=rules), backend=ManifestBackend(manifest), max_expansions=0)


@pytest.mark.parametrize("manifest", [
    "nonadjacent-pairs/model.json", "tethered-nac/model-guarded.json",
    "tethered-nac/model-ablated.json", "brusselator-growth/discovery-model.json",
    "local-slots/model-local.json", "gated-path/model.json",
])
def test_manifest_inputs_expose_actual_selected_raw_counts(manifest):
    from fractions import Fraction
    inputs = load(ROOT / "examples" / manifest)
    specification = inputs.specification()
    assert list(specification.observables) == inputs.data["initial_observables"]
    assert specification.omega.record()["definition"] == inputs.read("state_language")
    for value in specification.rules.values():
        specification.omega.validate_rule(value.rule)
    for query in specification.observables.values():
        coefficient, rule = query.identity_terms()[0]
        assert coefficient == Fraction(1)  # manifest normalization is readout only
        specification.omega.validate_rule(rule)


def test_manifest_analyze_checks_omega_before_worklist(monkeypatch):
    from mfobs.core import worklist
    inputs = load(ROOT / "examples/local-slots/model-local.json")
    omega = inputs.specification().omega
    calls = []
    monkeypatch.setattr(worklist, "run", lambda model, **options: calls.append(model) or {"status": "checked"})
    from mfobs import analyze, ManifestBackend
    backend = ManifestBackend(inputs)
    assert analyze(inputs.specification(), backend=backend, max_expansions=0).record() == {"status": "checked"}
    assert calls == [inputs]
    with pytest.raises(ValueError, match="unsupported Ω"):
        analyze(replace(inputs.specification(), omega=changed(omega, kind="different")), backend=backend)
    assert calls == [inputs]
