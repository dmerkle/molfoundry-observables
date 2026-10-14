from copy import deepcopy
import json

import pytest
import sympy as sp

from mfobs.invariants.compress import CompressionError, compress
from mfobs.outputs import matrix


def result(matrices, symbols=None, seeds=None):
    names = [f'F{i}' for i in range(next(iter(matrices.values())).rows)]
    return {
        'model': 'exact-linear-test', 'title': 'Exact linear test', 'closed': True,
        'basis': names, 'depths': {name: 0 if name in (seeds or [names[0]]) else 1 for name in names},
        'rates': [{'id': name, 'rate': (symbols or {}).get(name, name), 'per_match': '1', 'kind': 'abstract'}
                  for name in matrices],
        'equations_by_reaction': {
            source: {rate: {target: str(mat[i, j]) for j, target in enumerate(names) if mat[i, j]}
                     for rate, mat in matrices.items()}
            for i, source in enumerate(names)
        },
        'initial_values': {name: str(i+2) for i, name in enumerate(names)},
    }


def transform(compressed):
    return sp.Matrix(compressed['transformation_matrix'])


def test_exact_linear_combination_reduces_graph_representative_basis_and_preserves_source():
    a = sp.Matrix([[0, 1, 1, 0], [0, 0, 0, 1], [0, 0, 0, -1], [0, 0, 0, 0]])
    b = sp.Matrix([[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
    source = result({'a': a, 'b': b})
    saved = deepcopy(source)
    compressed = compress(source)
    assert source == saved
    assert compressed['status'] == 'EXACT_RATE_FAMILY_LINEAR_COMPRESSION'
    assert compressed['dimension'] == 2
    assert compressed['basis_combinations'] == {'K0': {'F0': '1'}, 'K1': {'F1': '1', 'F2': '1'}}
    assert compressed['initial_values'] == {'K0': '2', 'K1': '7'}
    assert compressed['seed_reconstruction'] == {'F0': {'K0': '1'}}
    assert transform(compressed) * matrix(source) == matrix(compressed) * transform(compressed)
    assert matrix(compressed) == sp.Matrix([[0, sp.Symbol('a')], [0, sp.Symbol('b')]])


def test_shared_rate_channels_are_aggregated_before_span_closure():
    plus = sp.Matrix([[0, 1], [0, 0]])
    source = result({'forward': plus, 'backward': -plus}, {'forward': 'k', 'backward': 'k'})
    # Stored coefficients already contain any per-match multipliers. Applying
    # this declaration a second time would destroy the exact cancellation.
    source['rates'][0]['per_match'] = '1/7'
    compressed = compress(source)
    assert compressed['dimension'] == 1
    assert compressed['aggregated_rate_channels'] == {'k': ['forward', 'backward']}
    assert matrix(compressed) == sp.zeros(1)
    # Distinct symbols remain independent, even if a later user chooses them
    # numerically equal. No generic parameter rank is guessed from sampling.
    independent = compress(result({'forward': plus, 'backward': -plus}))
    assert independent['dimension'] == 2


def test_multiple_explicit_seeds_exact_reconstruction_and_rational_coefficients():
    a = sp.Matrix([[0, sp.Rational(1, 3), 0], [0, 0, 0], [0, 0, 0]])
    source = result({'a': a})
    source['initial_observables'] = ['F2', 'F0']
    compressed = compress(source)
    assert compressed['dimension'] == 3
    assert compressed['initial_observables'] == ['K0', 'K1']
    assert compressed['source_initial_observables'] == ['F2', 'F0']
    t = transform(compressed)
    for seed, coefficients in compressed['seed_reconstruction'].items():
        reconstruction = sp.Matrix([[sp.Rational(coefficients.get(name, '0')) for name in compressed['basis']]])
        expected = sp.zeros(1, len(source['basis']))
        expected[0, source['basis'].index(seed)] = 1
        assert reconstruction * t == expected
    assert t * matrix(source) == matrix(compressed) * t
    assert compressed['basis_origins']['K2']['residual_scale'] == '1/3'


def test_rejects_unknown_incomplete_or_inexact_sources():
    source = result({'a': sp.zeros(2)})
    for patch in [{'closed': False}, {'unknown': [{'reason': 'budget'}]}, {'pending': ['F1']}, {'depths': {}}]:
        with pytest.raises(CompressionError):
            compress({**source, **patch})
    source['equations_by_reaction']['F0']['a']['F1'] = 0.5
    with pytest.raises(CompressionError, match='exact rational'):
        compress(source)
    source['equations_by_reaction']['F0']['a']['F1'] = 'k'
    with pytest.raises(CompressionError, match='exact rational'):
        compress(source)
    source['equations_by_reaction']['F0']['a']['F1'] = sp.Float('0.5')
    with pytest.raises(CompressionError, match='exact rational'):
        compress(source)


def test_absent_initial_values_not_assumed_zero_and_compression_is_repeatable():
    source = result({'a': sp.Matrix([[0, 1], [0, 0]])})
    del source['initial_values']
    first = compress(source)
    assert first['initial_values'] == {}
    assert first['initial_values_available'] is False
    second = compress(first)
    assert second['dimension'] == first['dimension']
    assert second['initial_values_available'] is False


def test_cli_exports_combinations_without_fabricating_graphs_and_preserves_source(tmp_path, capsys):
    from mfobs.cli import main
    source=result({'a':sp.Matrix([[0,1,1],[0,0,0],[0,0,0]])})
    path=tmp_path/'source.json';path.write_text(json.dumps(source), encoding='utf-8')
    original=path.read_bytes();out=tmp_path/'compressed'
    assert main(['compress',str(path),'--out',str(out)])==0
    compressed=json.loads((out/'result.json').read_text(encoding='utf-8'))
    assert compressed['dimension']==2
    assert (out/'transformation.csv').exists()
    assert json.loads((out/'seed-reconstruction.json').read_text(encoding='utf-8'))=={'F0':{'K0':'1'}}
    assert not (out/'identity-patterns').exists()
    assert 'rational linear combination' in (out/'ode.txt').read_text(encoding='utf-8')
    assert 'RAW injective pattern count' not in (out/'ode.txt').read_text(encoding='utf-8')
    assert path.read_bytes()==original
    assert main(['compress',str(path),'--out',str(out)])==1
    assert main(['compress',str(path),'--out',str(out),'--force'])==0
    # Even --force cannot overwrite the source result.
    assert main(['compress',str(out/'result.json'),'--out',str(out),'--force'])==1
    params=tmp_path/'parameters.json';params.write_text('{"a":"1/2"}', encoding='utf-8')
    assert main(['solve',str(out/'result.json'),'--parameters',str(params),
                 '--points','2','--until','1','--out',str(out/'trajectory.csv')])==0
    rows=(out/'trajectory.csv').read_text(encoding='utf-8').splitlines()
    assert rows[0]=='t,mu_K0,mu_K1'
    assert [float(v) for v in rows[-1].split(',')]==[1.0,5.5,7.0]
