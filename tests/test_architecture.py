import ast
from pathlib import Path

from mfobs.core.semantic_contract import SemanticContract
from mfobs.core.synthesis import DynamicBasisAdapter
from mfobs.inputs import load
from mfobs.models import eight_slots, chain_mixture, ring


ROOT=Path(__file__).resolve().parents[1]


def test_core_algorithm_is_model_independent():
    path=ROOT/'src/mfobs/core/linear_closure.py';text=path.read_text(encoding='utf-8')
    tree=ast.parse(text)
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom):assert not (node.module or '').startswith('mfobs.models')
        if isinstance(node,ast.Import):assert all(not item.name.startswith('mfobs.models') for item in node.names)


def test_composition_jump_compiler_and_incidence_operations_are_model_independent():
    for relative in (
        "src/mfobs/core/composition_jump_model.py",
        "src/mfobs/core/typed_incidence.py",
    ):
        tree = ast.parse((ROOT / relative).read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert not (node.module or "").startswith("mfobs.models")
            if isinstance(node, ast.Import):
                assert all(
                    not item.name.startswith("mfobs.models")
                    for item in node.names
                )


def test_initial_mixtures_are_component_lists_not_synthesis_bounds(tmp_path):
    from types import SimpleNamespace
    from mfobs.specifications import validate_initial_mixture
    components = []
    for label, count in [('A', 1), ('U', 1), ('T', 5)]:
        name = f'{label}.gml'
        (tmp_path / name).write_text(f'graph [ node [ id 0 label "{label}" ] ]', encoding="utf-8")
        components.append({'id': label, 'file': name, 'multiplicity': count})
    record = validate_initial_mixture(SimpleNamespace(file=tmp_path.joinpath), {
        'schema': 'explicit-component-mixture-v1',
        'description': 'Path, controller and spectator components.',
        'components': components,
    })
    assert [x['multiplicity'] for x in record['components']] == [1, 1, 5]


def test_productive_symbolic_models_share_the_executable_core():
    cases = [
        (
            chain_mixture.Adapter,
            'examples/gated-path/model.json',
            'theorem-argument',
        ),
        (
            eight_slots.Adapter,
            'examples/eight-slots/model.json',
            'machine-check',
        ),
    ]
    for adapter_type, manifest, evidence in cases:
        assert issubclass(adapter_type, DynamicBasisAdapter)
        inputs = load(ROOT / manifest)
        adapter = adapter_type(inputs)
        assert isinstance(adapter.semantic_contract, SemanticContract)
        contract = adapter.semantic_contract.record()
        assert contract['all_rules_verified']
        assert contract['preservation_evidence_levels'] == [evidence]
        assert set(contract['preservation']) == {
            rule['id'] for rule in inputs.data['rules']
        }


def test_fixed_ring_is_explicitly_an_exhaustive_reference_adapter():
    assert issubclass(ring.Adapter, DynamicBasisAdapter)
    # The control is shared; its reducer still enumerates the whole finite domain.
    assert 'exhaustive finite reference adapter' in ring.__doc__.lower()


def test_no_private_engine_access():
    allowed = set()
    found = set()
    for path in (ROOT / 'src/mfobs').rglob('*.py'):
        relative = path.relative_to(ROOT / 'src/mfobs').as_posix()
        tree = ast.parse(path.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and (node.module or '').startswith('molfoundry._'):
                found.add(relative)
            if isinstance(node, ast.Import) and any(item.name.startswith('molfoundry._') for item in node.names):
                found.add(relative)
        if relative not in allowed:
            assert '._ref' not in path.read_text(encoding='utf-8'), relative
    assert found == allowed
    for path in (ROOT / 'examples/tcs').glob('*.py'):
        assert '._ref' not in path.read_text(encoding='utf-8')
        assert 'molfoundry._engine' not in path.read_text(encoding='utf-8')


def test_restricted_composition_uses_only_public_molfoundry():
    for relative in ('core/composition.py', 'core/engine.py', 'core/forbidden.py',
                     'core/typed_incidence.py'):
        text = (ROOT / 'src/mfobs' / relative).read_text(encoding='utf-8')
        assert 'molfoundry._' not in text and '._ref' not in text, relative
    assert not (ROOT / 'src/mfobs/core/restricted_composition.py').exists()
    assert not (ROOT / 'src/mfobs/core/streaming.py').exists()
