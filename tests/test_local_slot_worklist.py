"""The local experiment must never query the finite Ring8 reduction oracle."""
from copy import deepcopy
from pathlib import Path

import pytest
import sympy as sp

from mfobs.core import mixtures
from mfobs.core.worklist import run
from mfobs.inputs import InputError, load
from mfobs.models import ring_language
from mfobs.outputs import export, matrix

ROOT = Path(__file__).resolve().parents[1]


def model(tag):
    return load(ROOT / f'examples/local-slots/model-{tag}.json')


@pytest.fixture
def no_host_oracle(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Synthesis queried a concrete-host/orbit evaluator')
    for name in ('orbit_generator', 'slot_ring_graph', 'evaluate_observable', 'word_values'):
        monkeypatch.setattr(ring_language, name, forbidden)
    monkeypatch.setattr(mixtures, 'fast_count', forbidden)


@pytest.mark.parametrize('tag', ['local', 'degree2'])
def test_first_step_comes_from_composition_without_oracle(tag, no_host_oracle):
    result = run(model(tag), max_expansions=1)
    assert not result['closed'] and result['matrix'] is None
    assert result['basis'] == ['P0', 'P1', 'P2']
    assert result['equations_by_reaction']['P0'] == {
        'basal_bind': {'P1': '1'}, 'cooperative_bind': {'P2': '1'},
        'unbind': {'P0': '-1'},
    }
    assert [result['normalized_patterns'][name]['quotient']['occupancies']
            for name in result['basis']] == [[1], [0], [0, 1]]
    assert [(r['left_raw'], r['right_raw']) for r in result['expansions']] == [
        (6, 8), (24, 27), (8, 6)]
    assert result['synthesis_state_evaluations'] == result['synthesis_orbit_evaluations'] == 0
    assert result['initial_values'] == {}


def test_noncooperative_system_closes_locally(no_host_oracle):
    result = run(model('no-coop'), max_expansions=3)
    assert result['closed'] and not result['pending'] and not result['unknown']
    assert result['status'] == 'CLOSED_MODULO_LOCAL_SLOT_AXIOMS'
    assert result['basis'] == ['P0', 'P1']
    k_on, k_off = sp.symbols('k_on k_off')
    assert matrix(result) == sp.Matrix([[-k_off, k_on], [k_off, -k_on]])
    assert all(proof['passed'] for proof in result['rule_preservation'].values())


@pytest.mark.parametrize('tag,size', [('local', 15), ('degree2', 13)])
def test_complete_first_three_waves_are_not_a_truncated_closed_ode(tag, size, tmp_path, no_host_oracle):
    result = run(model(tag), max_expansions=7)
    assert len(result['basis']) == size and len(result['steps']) == 7
    assert len(result['pending']) == size - 7
    assert result['unknown'] == [{'reason': 'expansion budget exhausted', 'pending': result['pending']}]
    assert result['status'] == 'UNKNOWN' and result['matrix'] is None
    export(result, tmp_path)
    assert (tmp_path / 'ode.txt').read_text(encoding='utf-8').startswith('UNKNOWN:')
    assert not (tmp_path / 'matrix.csv').exists()
    assert [r['after'] for r in result['rounds']] == [1, 3, 7, size]


@pytest.mark.parametrize('change', [{'connected': True}, {'include_empty': False}])
def test_restricted_overlap_shortcut_is_rejected(change):
    inputs = model('local')
    inputs.data = deepcopy(inputs.data)
    inputs.data['composition'].update(change)
    with pytest.raises(InputError):
        run(inputs, max_expansions=1)


def test_pattern_budget_is_unknown_not_a_semantic_zero(monkeypatch):
    inputs = model('local')
    original = inputs.read
    def read(key):
        value = original(key)
        return {**value, 'max_pattern_sites': 1} if key == 'reduction' else value
    monkeypatch.setattr(inputs, 'read', read)
    result = run(inputs, max_expansions=3)
    assert result['basis'] == ['P0', 'P1']
    assert not result['closed'] and result['matrix'] is None
    assert any('pattern site budget exceeded' in x['reason'] for x in result['unknown'])


def test_seed_name_cannot_collide_with_a_discovered_pattern(no_host_oracle):
    inputs = model('no-coop')
    inputs.data['initial_observables'] = ['P1']
    inputs.data['observables'][0]['id'] = 'P1'
    result = run(inputs, max_expansions=3)
    assert result['closed'] and result['basis'] == ['P1', 'P2']
    assert result['equations_by_reaction']['P1'] == {
        'basal_bind': {'P2': '1'}, 'unbind': {'P1': '-1'}}


@pytest.mark.parametrize('mode', ['public', 'restricted'])
def test_uncertified_external_pruning_cannot_produce_a_false_closure(mode):
    inputs = model('no-coop')
    inputs.data['forbidden'] = '../rop/forbidden.json'
    inputs.data['composition']['mode'] = mode
    with pytest.raises(InputError, match='pruning theories'):
        run(inputs, max_expansions=3)
