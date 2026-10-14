from mfobs.core.graph import make_graph, edge_rows, matches as plain_matches
from mfobs.core.rules import make_rule
from pathlib import Path
import json,ast
import pytest
from mfobs.inputs import load,InputError
from mfobs.core.conditions import Adjacency, matches, UnsupportedCondition
from mfobs.core.worklist import run
from mfobs.outputs import export,ode_text
from mfobs.core.engine import mf
ROOT=Path(__file__).resolve().parents[1]

def make_birth(tmp_path,mode='public',factor='1'):
    files={'birth.gml':'rule [ ruleID "birth" right [ node [ id 0 label "A" ] ] ]',
           'count.gml':'rule [ ruleID "count-A" context [ node [ id 0 label "A" ] ] ]'}
    for name,text in files.items():(tmp_path/name).write_text(text, encoding='utf-8')
    for name,data in [('state.json',{'scope':'all finite A-labelled graphs'}),('reduction.json',{'method':'canonical only'}),('parameters.json',{'beta':'1'})]:
        (tmp_path/name).write_text(json.dumps(data), encoding='utf-8')
    d={'schema_version':1,'id':'birth_test','title':'Independent test fixture','adapter':'mfobs.core.identity:Adapter',
       'rules':[{'id':'birth','file':'birth.gml','rate':'beta','per_match':factor,'kind':'abstract'}],
       'observables':[{'id':'N','file':'count.gml'}],'initial_observables':['N'],
       'composition':{'mode':mode,'connected':False,'include_empty':True},'forbidden':None,
       'state_language':'state.json','reduction':'reduction.json','normalisation':{},'initial_state':{},
       'parameters':'parameters.json','reference':None}
    p=tmp_path/'model.json';p.write_text(json.dumps(d), encoding='utf-8');return load(p)

def test_generic_birth_closes_without_model_reducer(tmp_path):
    result=run(make_birth(tmp_path))
    assert result['closed'] and len(result['basis'])==2
    assert result['equations_by_reaction']['N']['birth']=={result['basis'][1]:'1'}
    assert not result['equations_by_reaction'][result['basis'][1]]['birth']


def test_raw_per_match_multiplier_is_input(tmp_path):
    result=run(make_birth(tmp_path,factor='2'))
    assert result['equations_by_reaction']['N']['birth']=={result['basis'][1]:'2'}


def test_expansion_budget_is_unknown_not_false_nonclosure(tmp_path):
    checkpoint = tmp_path/'checkpoint.json'
    result=run(make_birth(tmp_path),max_expansions=1,checkpoint=checkpoint)
    assert not result['closed'] and result['status']=='UNKNOWN' and result['pending']
    assert result['matrix'] is None
    assert ode_text(result).startswith('UNKNOWN')
    partial = json.loads(checkpoint.read_text(encoding='utf-8'))
    assert partial['basis'] == result['basis']
    assert partial['pending'] == result['pending']


def test_unknown_overwrite_removes_stale_ode(tmp_path):
    inp=make_birth(tmp_path);out=tmp_path/'out';export(run(inp),out)
    assert (out/'ode.tex').exists()
    import csv
    with (out/'coefficients-by-rate.csv').open(encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]['rate_channel'] == 'beta'
    assert rows[0]['coefficient_after_per_match_factor'] == '1'
    (out/'verification.json').write_text('{"passed":true}', encoding='utf-8')
    (out/'trajectory.csv').write_text('old trajectory', encoding='utf-8')
    export(run(inp,max_expansions=0),out)
    assert not (out/'ode.tex').exists() and not (out/'matrix.csv').exists()
    assert not (out/'verification.json').exists() and not (out/'trajectory.csv').exists()

@pytest.mark.parametrize('model',['gated-path','ring8'])
def test_reviewed_inputs_load_and_match_contract(model):
    inp=load(ROOT/'examples'/model/'model.json');inp.adapter()


def test_invalid_manifest_fields_fail(tmp_path):
    inp=make_birth(tmp_path);d=inp.data;d['silently_ignored_field']=1;inp.path.write_text(json.dumps(d), encoding='utf-8')
    with pytest.raises(InputError):load(inp.path)


def test_duplicate_rules_fail(tmp_path):
    inp=make_birth(tmp_path);d=inp.data;d['rules']*=2;inp.path.write_text(json.dumps(d), encoding='utf-8')
    with pytest.raises(InputError):load(inp.path)


def test_core_has_no_model_imports():
    for p in (ROOT/'src/mfobs/core').glob('*.py'):
        for n in ast.walk(ast.parse(p.read_text(encoding='utf-8'))):
            if isinstance(n,ast.ImportFrom):assert not (n.module or '').startswith('mfobs.models'),p
            if isinstance(n,ast.Import):assert all(not a.name.startswith('mfobs.models') for a in n.names),p


def test_full_host_neighbors_not_only_matched_neighbors():
    pattern=make_graph({0:'C',1:'S'},[(0,1,'-')])
    host=make_graph({0:'C',1:'S',2:'C',3:'C',4:'C'},[(0,1,'-'),(0,2,'-'),(0,3,'-'),(0,4,'-')])
    assert len(plain_matches(pattern,host))==1
    assert not matches(pattern,host,adjacency=(Adjacency(0,('C',),('-',),'<',3),))

@pytest.mark.parametrize('op',['!=','=='])
def test_unknown_operators_rejected(op):
    with pytest.raises(UnsupportedCondition):Adjacency(0,('C',),('-',),op,1)


def test_unsupported_condition_returns_unknown():
    rule=mf.Rule.fromGMLString('rule [ context [ node [id 0 label "A"] node [id 1 label "A"] ] constrainShortestPath [ source 0 target 1 op "<" length 3 ] ]',add=False)
    from mfobs.invariants.counting import Reducer
    reducer=Reducer(load(ROOT/'examples/gated-path/model.json').read('state_language'))
    assert reducer.reduce(rule.canonKey,rule)['status']=='unknown'


def test_partial_pushout_invalidity_not_monotone_in_overlap():
    partial=make_graph({0:'H',1:'C',2:'C'},[(0,1,'-'),(0,2,'-')])
    full=make_graph({0:'H',1:'C'},[(0,1,'-')])
    degree=lambda g:sum(0 in (u,v) for u,v,_ in edge_rows(g))
    assert degree(partial)==2 and degree(full)==1


@pytest.mark.parametrize('value',[0.5,'0.5','-1','0','1/0'])
def test_noncanonical_rate_factor_rejected(tmp_path,value):
    inp=make_birth(tmp_path);inp.data['rules'][0]['per_match']=value
    inp.path.write_text(json.dumps(inp.data), encoding='utf-8')
    with pytest.raises(InputError):load(inp.path)

def test_adapter_error_does_not_change_later_compositions(tmp_path):
    inp=make_birth(tmp_path);inp.data['adapter']='does_not_exist:Adapter'
    with pytest.raises(ModuleNotFoundError):run(inp)
    assert run(make_birth(tmp_path))['closed']


def test_numerical_rate_values_not_read_by_synthesis(tmp_path):
    inp=make_birth(tmp_path)
    (tmp_path/'parameters.json').write_text('not numeric JSON; only solve should care', encoding='utf-8')
    result=run(load(inp.path))
    assert result['closed'] and result['equations_by_reaction']['N']['birth']=={'O1':'1'}


def test_injective_embedding_multiplicity_not_deduplicated():
    from mfobs.core.graph import observable
    pattern=make_graph({0:'A',1:'A'})
    host=make_graph({0:'A',1:'A',2:'A'})
    assert observable(pattern,host)==6


def test_dpo_dangling_deletion_is_not_spo_deletion():
    from mfobs.core.graph import transitions
    rule=make_rule('delete-A',make_graph({0:'A'}),make_graph({}))
    host=make_graph({0:'A',1:'B'},[(0,1,'e')])
    assert not transitions(rule,host,'DPO')
    with pytest.raises(ValueError, match='Only MolFoundry DPO'):
        transitions(rule,host,'SPO')


def test_relational_no_edge_has_actual_graph_semantics():
    from mfobs.core.graph import transitions
    left=make_graph({0:'A',1:'B'})
    right=make_graph({0:'A',1:'B'},[(0,1,'x')])
    rule=make_rule('bind',left,right,((0,1),))
    assert sum(transitions(rule,left,'DPO').values())==1
    assert not transitions(rule,right,'DPO')


def test_cli_available_from_outside_checkout(tmp_path):
    import subprocess,sys
    proc=subprocess.run([sys.executable,str(ROOT/'run.py'),'inspect',str(ROOT/'examples/ring8/model.json')],cwd=tmp_path,capture_output=True,text=True)
    assert proc.returncode==0,proc.stderr
    assert 'k_on' in proc.stdout and 'Initially queued: P0' in proc.stdout
