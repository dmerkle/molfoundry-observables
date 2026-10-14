"""Scientific checks of fixed-ring walkthrough data."""
from fractions import Fraction
import gzip
import json
import os
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as sp

from mfobs.inputs import load
from mfobs.models import ring_language as ring
from mfobs.outputs import matrix
from conftest import regression_record

ROOT = Path(os.environ.get('MF_WALKTHROUGH_ROOT', Path(__file__).resolve().parents[1]))
RING_DATA = Path(os.environ.get('MF_RING_WALKTHROUGH_DATA', ROOT/'experiments/ring8/walkthrough-data.json'))


def read(path):
    path = Path(path)
    if path.suffix == '.gz':
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            return json.load(stream)
    return json.loads(path.read_text(encoding='utf-8'))


def recorded_path(path):
    assert not Path(path).is_absolute(), f'Unexpected nonportable source: {path}'
    return ROOT/path


@pytest.fixture(scope='module')
def ring_data():
    return regression_record('ring-dynamics')


def test_ring_provenance_and_exact_normalized_rate_matrices(ring_data):
    d=ring_data; result=regression_record('ring-result')
    assert d['input_definition_sha256']==load(ROOT/'examples/ring8/model.json').definition_digest()==result['input_definition_sha256']
    assert result['closed'] and not result['unknown'] and len(d['basis'])==30
    assert d['computed_basis']==result['basis'] and d['scales']==[8,8]+[16]*28
    assert d['basis']==[f'M{i}' for i in range(30)]
    assert len(d['orbits'])==30 and len(d['patterns'])==30 and len(d['rules'])==3
    assert d['discovery']['rounds']==result['rounds']
    assert d['discovery']['raw_compositions']==result['raw_compositions']
    assert d['discovery']['ordered_composition_calls']==result['ordered_composition_calls']
    b=sp.Matrix([[result['basis_vectors'][name][i] for i in range(30)] for name in d['computed_basis']])
    transform=sp.Matrix(d['presentation_from_computed'])
    assert transform.det()!=0
    b=transform*b
    assert b==sp.Matrix([p['raw_values_on_orbits'] for p in d['patterns']])
    scale=sp.diag(*d['scales']); raw=transform*matrix(result)*transform.inv()
    for channel, symbol in d['matrices']['rate_symbols'].items():
        a=raw.subs({sp.Symbol(rate):int(rate==symbol) for rate in d['matrices']['rate_symbols'].values()})
        assert sp.Matrix(d['matrices']['raw_by_rate'][channel])==a
        assert sp.Matrix(d['matrices']['normalized_by_rate'][channel])==scale.inv()*a*scale
        _,q=ring.orbit_generator(8,ring.RingRates(**{k:int(k==channel) for k in d['matrices']['rate_symbols']}))
        assert a*b==b*q.T
    v=d['verification']
    assert v['passed'] and v['independent_graph_pattern_counts']==900
    assert v['independent_bit_motif_counts']==7680 and v['rate_separated_full_state_rows']==23040
    assert v['saved_multiplicity_comparisons']==regression_record('ring-verification')['multiplicities']['comparisons']
    assert v['maximum_normalized_ODE_reference_error']<1e-9


def test_ring_initial_vectors_spatial_rate_comparison_and_sampled_solutions(ring_data):
    d=ring_data; cases={x['id']:x for x in d['cases']}; scenarios={x['id']:x for x in d['scenarios']}
    assert [cases[k]['bits'] for k in ('empty','block4','alternating4')]==[[0]*8,[1,1,1,1,0,0,0,0],[1,0]*4]
    assert {k:cases[k]['raw_channel_matches'] for k in cases}=={
        'empty':{'basal_bind':8,'cooperative_bind':0,'unbind':0},
        'block4':{'basal_bind':4,'cooperative_bind':2,'unbind':4},
        'alternating4':{'basal_bind':4,'cooperative_bind':8,'unbind':4}}
    for case in cases.values():
        assert len(case['initial_normalized'])==len(case['initial_raw'])==30
        assert [Fraction(v)/s for v,s in zip(case['initial_raw'],d['scales'])]==list(map(Fraction,case['initial_normalized']))
        assert len(case['graph']['nodes'])==24
    assert len(d['curves'])==9
    curves={(x['scenario'],x['case']):x for x in d['curves']}
    rate_matrices={channel:np.array([list(map(Fraction,row)) for row in rows],float)
                   for channel,rows in d['matrices']['normalized_by_rate'].items()}
    for (scenario,initial),curve in curves.items():
        rates=scenarios[scenario]['parameters'];a=np.zeros((30,30))
        for channel,symbol in d['matrices']['rate_symbols'].items():
            a+=float(Fraction(rates[symbol]))*rate_matrices[channel]
        y0=np.array(list(map(Fraction,cases[initial]['initial_normalized'])),float)
        assert len(curve['rows'])==241 and curve['rows'][0]['t']==0 and curve['rows'][-1]['t']==12
        for row in curve['rows'][::60]:
            np.testing.assert_allclose(row['values'],expm(row['t']*a)@y0,atol=1e-10,rtol=1e-10)
    assert Fraction(curves['cooperative','block4']['initial_derivative'][0])*8==Fraction(8,5)
    assert Fraction(curves['cooperative','alternating4']['initial_derivative'][0])*8==Fraction(44,5)
    for left,right in zip(curves['independent','block4']['rows'],curves['independent','alternating4']['rows']):
        assert abs(left['values'][0]-right['values'][0])<1e-10
    assert all(row['values'][0]==0 for row in curves['washout','empty']['rows'])


def test_ring_direct_expansion_is_bounded_graph_application_not_composition(ring_data):
    d=ring_data['direct_expansion'];graphs={g['id']:g for g in d['graphs']}
    assert len(graphs)==28 and len({g['canonical_key'] for g in graphs.values()})==28
    assert d['calls']==d['independent_match_checks']==d['independent_successor_set_checks']==93
    assert {c['id']:(c['state_count'],c['edge_count']) for c in d['cases']}=={
        'empty':(11,29),'block4':(24,116),'alternating4':(19,61)}
    for c in d['cases']:
        assert c['depth_limit']==3 and len(c['states'])==c['state_count'] and len(c['edges'])==c['edge_count']
        ids={s['graph'] for s in c['states']};assert len(ids)==len(c['states'])
        assert all(0<=s['depth']<=3 for s in c['states'])
        for edge in c['edges']:
            assert edge['source'] in ids and edge['target'] in ids and edge['raw_match_multiplicity']>0
            change=graphs[edge['target']]['occupied']-graphs[edge['source']]['occupied']
            assert change==(-1 if edge['rule']=='unbind' else 1)


def test_ring_reference_rejects_changed_generator():
    from mfobs.checks.reference import verify
    inputs=load(ROOT/'examples/ring8/model.json')
    result=regression_record('ring-result')
    assert verify(inputs,result)['passed']
    # A coordinate change is allowed; a changed rate-channel coefficient is not.
    row=result['equations_by_channel']['P0']['k_on']
    row['P0']=str(Fraction(row.get('P0','0'))+1)
    with pytest.raises(AssertionError,match='generator changed'):
        verify(inputs,result)
