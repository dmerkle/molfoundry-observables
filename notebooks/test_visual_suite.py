
import mfobs.models.remy.analysis as remy_analysis
import mfobs.models.dna_walker.inputs as walker_model_inputs
from pathlib import Path
import json,sys
import numpy as np
import nbformat
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT),str(ROOT.parent/'src'),str(ROOT.parent)]
import visual_suite as suite
from mfobs.checks import dna_walker as walker
from mfobs.core.graph import matches
from mfobs.core.rules import rule_graph


def test_recovered_walker_interfaces_and_parallel_match_multiplicity():
    evidence=ROOT/'visual/evidence/dna-walker'
    record=json.loads((evidence/'result.json').read_text(encoding='utf-8'))
    provenance=json.loads((evidence/'provenance.json').read_text(encoding='utf-8'))
    original,incidence=walker.verified_exports(evidence/'recovered',record,provenance)
    for rule in original['rules']:
        converted=walker_model_inputs.convert_rule(rule)
        # One leg edge is replaced, the other and the DNA edge are preserved.
        assert sum(v.leftLabel is not None for v in converted.vertices)==6
        assert sum(v.rightLabel is not None for v in converted.vertices)==6
        assert sum(v.membership.value=='context' for v in converted.vertices)==5
        assert not converted.matchConstraintsNoEdge
    compressed=rule_graph(walker_model_inputs.identity(original['compressed_query'],'Craw'))
    extended=rule_graph(walker_model_inputs.identity(original['extended_query'],'E'))
    assert len(matches(compressed,compressed))==2
    assert len(matches(extended,extended))==1
    assert all(row['Equal'] for row in walker.compare_incidence(record,incidence))


def test_walker_generated_probability_matrix_and_stationary_state():
    r=json.loads((ROOT/'visual/evidence/dna-walker/result.json').read_text(encoding='utf-8'))
    rates=dict(kFE=2,kBE=.4,kFC=1.5,kBC=.3)
    names,A=suite.walker_matrix(r,rates)
    stationary=np.array([2.4,1.8])/4.2
    assert names==['E','C']
    assert np.allclose(A@stationary,0)
    assert np.allclose(A.sum(axis=0),0)


def test_rate_controls_keep_horizon_runs_seed_and_solver_settings(monkeypatch):
    displays=[];calls=[]
    monkeypatch.setattr(suite,'display',displays.append)
    opts=dict(end=100.,paths=256,seed=72,truncation_tolerance=1e-8,max_states=12345)
    state=suite.rate_controls(lambda rates,**options:calls.append((rates,options)),{'rate':.5},opts)
    box=displays[0];box.children[0].value=.9;box.children[-2].click()
    assert calls==[({'rate':.9},opts)]
    assert state['options']==opts


def test_notebooks_use_portable_kernel_and_matching_sources():
    import re
    for p in (ROOT/'examples').glob('*.ipynb'):
        if not p.name[:2].isdigit() or int(p.name[:2])<2:continue
        n=nbformat.read(p,as_version=4);nbformat.validate(n)
        assert n.metadata.kernelspec.name=='molfoundry'
        script=p.with_suffix('.py').read_text(encoding='utf-8')
        blocks=re.findall(r'^# %%\n(.*?)(?=^# %%|\Z)',script,flags=re.M | re.S)
        cells=[c for c in n.cells if c.cell_type=='code']
        assert len(cells)==len(blocks)
        for c,block in zip(cells,blocks):
            compile(c.source,str(p),'exec')
            assert c.source.strip()==block.strip()


def test_remy_exported_gml_is_the_discovered_expression():
    from hashlib import sha256
    import molfoundry as mf
    record=remy_analysis.common_analysis(order=3,linear_limit=16,record_patterns=True)
    for expression in record['linear']['raw_rule_expressions'].values():
        assert len(expression['gml_terms'])==expression['term_count']
        computed=[]
        for term in expression['gml_terms']:
            rule=mf.Rule.fromGMLString(term['gml'],add=False)
            computed.append({'term_id':sha256(rule.canonKey.encode()).hexdigest(),
                             'coefficient':term['coefficient']})
        assert computed==expression['terms']


def test_examples_do_not_run_post_synthesis_checks():
    for p in [*(ROOT/'examples').glob('*.ipynb'),
              *(ROOT/'visual/comparisons').glob('*.ipynb')]:
        book = nbformat.read(p, as_version=4)
        text = '\n'.join(c.source for c in book.cells)
        assert 'visual.check(' not in text, p
        assert 'RUN_CHECKS' not in text, p
        assert 'mfobs.checks' not in text, p
        assert 'id="checks"' not in text, p


def test_brusselator_approximation_conserves_monomer_plus_chain_balance():
    from scipy.integrate import solve_ivp
    from mfobs.models.brusselator_growth_chemistry import deterministic_rhs
    initial = [.4, 2.5, 4., 0.]
    rates = dict(a=1., b=2.2, k1=1., k2=1., k3=1., k4=1., beta=.3, k_g=.15)
    times = np.linspace(0, 60, 101)
    solution = solve_ivp(lambda t, y: deterministic_rhs(t, y, rates),
                         (0, 60), initial, t_eval=times, rtol=1e-9, atol=1e-11)
    assert solution.success
    assert np.allclose(solution.y[2]+solution.y[3], initial[2]+initial[3]+rates['beta']*times, atol=1e-7)


def test_growing_ring_approximation_conserves_nonnegative_site_fractions():
    from scipy.integrate import solve_ivp
    from mfobs.models.structural_nac_oscillator import one_site_rhs
    for growth in (0., .12):
        solution = solve_ivp(lambda t, y: one_site_rhs(y, 1., growth),
                             (0, 80), [.72, .14, .14], rtol=1e-9, atol=1e-11)
        assert solution.success
        assert np.allclose(solution.y.sum(axis=0), 1, atol=1e-9)
        assert solution.y.min() >= -1e-9
