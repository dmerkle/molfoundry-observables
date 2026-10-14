from mfobs.core.graph import node_rows, make_graph, edge_rows, node_labels
from pathlib import Path
from fractions import Fraction
import ast,json,random,shutil
import networkx as nx
import pytest

from mfobs.core.conditions import Adjacency
from mfobs.core.mixtures import fast_count
from mfobs.core.gml import disjoint_union
from mfobs.core.engine import mf
from mfobs.invariants.chain import ChainLanguage,ChainCounter
from mfobs.invariants.sequences import Sequence,footprint_sequence
from mfobs.invariants.counting import Reducer,evaluate
from mfobs.core.linear_closure import RationalSpan as Span
from mfobs.inputs import load
from mfobs.core.worklist import run
ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def language():return json.loads((ROOT/'examples/gated-path/state-language.json').read_text(encoding='utf-8'))

@pytest.fixture(scope='module')
def exact(language):return Reducer(language)

def nxg(g):
    h=nx.Graph();h.add_nodes_from((u,{'l':l}) for u,l in node_rows(g));h.add_edges_from((u,v,{'b':b}) for u,v,b in edge_rows(g));return h

def test_independent_grammar_matches_paths(language):
    lang=ChainLanguage(language['chain'])
    for n in [0,1,2,3,5,10,25,40]:
        path = make_graph({i: 'A' for i in range(n+1)}, [(i, i+1, '-') for i in range(n)])
        a,b=nxg(path),nxg(lang.instantiate(n)[0])
        assert nx.is_isomorphic(a,b,node_match=lambda a,b:a==b,edge_match=lambda a,b:a==b)

@pytest.mark.parametrize('length',[1,2,3,5,8])
def test_footprint_theorem_on_ordinary_path_graphs(length):
    # Capped one-vertex tiles give P_(n+1). The count is known independently.
    g={'nodes':[[0,'C']],'edges':[]}
    lang=ChainLanguage({'kind':'one-port-decorated-path','join_bond':'-', 'zero':g,
        'left':{**g,'incoming':None,'outgoing':0},'repeat':{**g,'incoming':0,'outgoing':0},'right':{**g,'incoming':0,'outgoing':None}})
    p=make_graph({i:'C' for i in range(length)},[(i,i+1,'-') for i in range(length-1)])
    s=ChainCounter(lang).count(p).sequence
    for n in [*range(15),1000]:
        want=n+1 if length==1 else 2*max(0,n+2-length)
        assert s.value(n)==want


def test_random_connected_footprints_and_conditions(exact):
    rng=random.Random(91872);language=exact.chain.language
    for trial in range(60):
        host,_=language.instantiate(rng.randrange(5));adj={u:set() for u in node_labels(host)}
        for u,v,b in edge_rows(host):adj[u].add(v);adj[v].add(u)
        ids={rng.choice(list(node_labels(host)))};desired=rng.randrange(1,9)
        while len(ids)<desired:
            boundary=set().union(*(adj[u] for u in ids))-ids
            if not boundary:break
            ids.add(rng.choice(sorted(boundary)))
        p=make_graph({u:l for u,l in node_rows(host) if u in ids},
                     [(u,v,b) for u,v,b in edge_rows(host) if u in ids and v in ids]);ne=[];conditions=[]
        if len(ids)>=2 and trial%3==0:
            u,v=rng.sample(sorted(ids),2);ne=[(u,v)] # Includes contradictions.
        if trial%4==0:conditions=[Adjacency(min(ids),('A',),('-',),'<',rng.randrange(1,5))]
        result=exact.chain.count(p,tuple(ne),tuple(conditions))
        for n in [0,1,2,4,7,12]:
            assert result.sequence.value(n)==fast_count(p,language.instantiate(n)[0],ne,conditions),(trial,n,result.witness)


def test_disconnected_joint_count_is_exact(exact):
    p=make_graph({0:'A',1:'A'})
    result=exact.count(p)
    assert result['status']=='exact'
    for n,m,x in [(0,0,0),(1,2,1),(5,3,0)]:
        total = n + 1
        assert evaluate(result['polynomial'],n,m,x)==total*(total-1)


@pytest.mark.parametrize('copies',[1,2,3,4])
def test_falling_factorials_are_not_limited_to_degree_one(exact,copies):
    p=make_graph({0: 'T'})
    pp=disjoint_union([p]*copies);result=exact.count(pp)
    assert result['status']=='exact' and result['polynomial']=={('power',0,copies,0):Fraction(1)}
    for m in range(8):
        ff=1
        for j in range(copies):ff*=m-j
        assert evaluate(result['polynomial'],333,m,1)==ff


def test_derived_atom_caps(exact):
    assert exact.caps=={'U':1,'V':1,'W':1}
    p=make_graph({0:'U',1:'U'});r=exact.count(p)
    assert r['status']=='zero' and r['witness']['method']=='derived-atom-cap'


def test_unhandled_condition_is_unknown(exact):
    rule=mf.Rule.fromGMLString('rule [ context [ node [id 0 label "A"] node [id 1 label "A"] ] constrainShortestPath [ source 0 target 1 op "<" length 3 ] ]',add=False)
    assert exact.reduce(rule.canonKey,rule)['status']=='unknown'


def test_invariants_do_not_import_model_equations():
    paths=list((ROOT/'src/mfobs/invariants').glob('*.py'))+[ROOT/'src/mfobs/models/chain_mixture.py']
    for path in paths:
        for n in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            if isinstance(n,ast.ImportFrom):assert not (n.module or '').startswith('mfobs.models'),path
            if isinstance(n,ast.Import):assert all(not a.name.startswith('mfobs.models') for a in n.names),path


def test_no_prespecified_dimension_in_linear_algebra():
    span=Span()
    for q in range(12):assert span.insert(f'F{q}',{('power',0,q,0):Fraction(1)})
    residual,coords=span.reduce({('power',0,11,0):Fraction(3),('power',0,4,0):Fraction(-2)})
    assert not residual and coords=={'F4':-2,'F11':3}


def test_discovery_budget_does_not_claim_closure():
    result=run(load(ROOT/'examples/gated-path/model.json'),max_expansions=1)
    assert not result['closed'] and result['status']=='UNKNOWN' and result['matrix'] is None


def test_sequence_product_keeps_point_corrections():
    a=Sequence.make([1,2],{0:3,2:-5});b=Sequence.make([0,1],{1:2,2:7})
    for n in range(10):assert (a*b).value(n)==a.value(n)*b.value(n)


def test_footprint_counts_replayed_from_classes(exact):
    result=exact.chain.count(make_graph({0:'A',1:'A'}, [(0,1,'-')]))
    assert footprint_sequence(result.witness)==Sequence.make([0,2])


def test_discovery_without_reference_or_supplied_target_basis(tmp_path):
    destination=tmp_path/'example'
    shutil.copytree(ROOT/'examples/gated-path',destination)
    path=destination/'model.json'
    manifest=json.loads(path.read_text(encoding='utf-8'))
    manifest['reference']='missing-reference.json'
    path.write_text(json.dumps(manifest), encoding='utf-8')
    result=run(load(path))
    assert result['closed'] and len(result['basis'])==3 and not result['unknown']


def test_nac_ablation_changes_discovered_dimension():
    guarded=run(load(ROOT/'examples/gated-path/model.json'))
    ablated_inputs=load(ROOT/'examples/gated-path/model-ungated.json')
    ablated=run(ablated_inputs)
    assert guarded['closed'] and ablated['closed']
    assert len(guarded['basis'])==3 and len(ablated['basis'])==2
    from mfobs.checks.invariant_discovery import verify
    assert verify(ablated_inputs,ablated)['passed']
