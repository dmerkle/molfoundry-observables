#!/usr/bin/env python3
"""Independent finite-host verification of completed local-slot generator rows.

This checker never calls synthesis or its normalizer. Actual emitted GML is
parsed independently and counted by a generic typed injective matcher. These
finite tests are distinct from the accompanying all-cycle-size argument.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import re
import sys
from time import perf_counter

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mfobs.checks import local_slots_reference as ref
from mfobs.inputs import load


def parse_pattern(gml):
    """Strict independent parser for supported vertex-preserving jump identities."""
    tokens=re.findall(r'"(?:\\.|[^"\\])*"|\[|\]|[^\s\[\]]+',gml);position=0
    def block(nested=False):
        nonlocal position
        result={}
        while position<len(tokens) and tokens[position]!=']':
            key=tokens[position];position+=1
            if position==len(tokens):raise ValueError('Missing GML value')
            token=tokens[position];position+=1
            value=block(True) if token=='[' else json.loads(token) if token.startswith('"') else int(token)
            result.setdefault(key,[]).append(value)
        if nested:
            if position==len(tokens):raise ValueError('Unclosed GML block')
            position+=1
        return result
    document=block()
    assert position==len(tokens) and set(document)=={'rule'} and len(document['rule'])==1
    rule=document['rule'][0]
    assert not set(rule)-{'ruleID','context','left','right','constrainNoEdge'}
    for side in ('left','right'):
        assert not any(piece.get('node') for piece in rule.get(side,[])), 'Vertex-deleting jumps need a gluing checker'
    def side(which):
        pieces=rule.get('context',[])+rule.get(which,[])
        nodes=tuple(sorted((n['id'][0],n['label'][0]) for p in pieces for n in p.get('node',[])))
        edges=tuple(sorted((min(e['source'][0],e['target'][0]),max(e['source'][0],e['target'][0]),e['label'][0])
                           for p in pieces for e in p.get('edge',[])))
        assert len(dict(nodes))==len(nodes)
        assert all(u in dict(nodes) and v in dict(nodes) and u!=v for u,v,_ in edges)
        return ref.TypedGraph(nodes,edges)
    graph=side('left');assert graph==side('right'), 'Only diagonal jump patterns are counted'
    no_edges=tuple((c['source'][0],c['target'][0]) for c in rule.get('constrainNoEdge',[]))
    assert all(u in dict(graph.nodes) and v in dict(graph.nodes) for u,v in no_edges)
    return graph,no_edges


def host_families(max_degree):
    families=[('empty',0,()),('single',1,()),('discrete4',4,()),
              ('path4',4,((0,1),(1,2),(2,3))),('edge-and-isolate',3,((0,1),)),
              ('triangle-and-isolate',4,ref.cycle(3)),
              *[(f'C{n}',n,ref.cycle(n)) for n in (3,4,5,6,8)],
              ('star4',4,((0,1),(0,2),(0,3))),('K4',4,tuple(combinations(range(4),2)))]
    for name,n,edges in families:
        if max_degree is not None and max(map(len,ref.neighbours(n,edges)),default=0)>max_degree:
            continue
        # All assignments through six sites. Eight-site hosts use a bounded
        # deterministic sample plus all one-step successors for exact Qf.
        masks=list(range(1<<n)) if n<=6 else [0,255,15,85,1,3,7,129]
        states={tuple((mask>>i)&1 for i in range(n)) for mask in masks}
        closure=set(states)
        for bits in states:
            closure.update(target for target,_,_ in ref.events(bits,edges))
        yield name,n,edges,sorted(states),sorted(closure)


def verify_variant(variant):
    path=ROOT/f'out/local-slots/{variant}/result.json'
    result=json.loads(path.read_text(encoding='utf-8'));inputs=load(ROOT/f'examples/local-slots/model-{variant}.json')
    assert result['input_definition_sha256']==inputs.definition_digest()
    assert result['synthesis_state_evaluations']==result['synthesis_orbit_evaluations']==0
    for source,expected in result['source_sha256'].items():
        assert sha256((ROOT/source).read_bytes()).hexdigest()==expected
    rates={row['id']:row for row in result['rates']}
    assert all(Fraction(r['per_match'])==1 for r in rates.values())
    channels={'basal_bind':'basal','cooperative_bind':'cooperative','unbind':'off'}
    parsed={name:parse_pattern(record['gml']) for name,record in result['identity_patterns'].items()}
    excluded={r['source'] for r in result['unknown'] if 'source' in r}
    complete=[step['expanded'] for step in result['steps'] if step['expanded'] not in excluded]
    assert complete and all(set(result['equations_by_reaction'][name])==set(rates) for name in complete)
    quotient_checks=0;generator_checks=0;families=[]
    for family,n,edges,states,closure in host_families(result['local_contract']['next_max_degree']):
        counts={}
        for bits in closure:
            host=ref.typed_slots(bits,edges)
            values={name:ref.typed_match_count(graph,host,ne) for name,(graph,ne) in parsed.items()}
            for name,value in values.items():
                quotient=result['normalized_patterns'][name]['quotient']
                assert value==ref.motif_count(bits,edges,quotient['occupancies'],quotient['next_edges'],quotient['no_next_edges']), (variant,family,bits,name)
                quotient_checks+=1
            counts[bits]=values
        records=[]
        for bits in states:
            for name in complete:
                for channel in rates:
                    params={key:int(key==channels[channel]) for key in ('basal','cooperative','off')}
                    direct=ref.generator_value(lambda state:counts[state][name],bits,edges,**params)
                    row=result['equations_by_reaction'][name][channel]
                    predicted=sum(Fraction(c)*counts[bits][other] for other,c in row.items())
                    assert direct==predicted,(variant,family,bits,name,channel,direct,predicted)
                    generator_checks+=1
                    if name=='P0':
                        records.append({'bits':list(bits),'channel':channel,'derivative':str(direct)})
        families.append({'name':family,'sites':n,'next_edges':[list(e) for e in edges],
            'evaluated_source_assignments':len(states),'counted_assignments_including_successors':len(closure),
            'occupancy_assignments_exhaustive_for_this_graph':len(states)==2**n,
            'first_equation_checks':records})
    status=result['status']
    return {'passed':True,'scope':'FINITE_HOST_VERIFICATION_ONLY',
        'result_status':status,'result_closed':result['closed'],'pending':result['pending'],
        'source_result':path.relative_to(ROOT).as_posix(),'source_result_sha256':sha256(path.read_bytes()).hexdigest(),
        'input_definition_sha256':inputs.definition_digest(),
        'patterns_counted':list(parsed),'completed_source_patterns_checked':complete,
        'unsupported_source_patterns_not_checked':sorted(excluded),
        'independent_typed_vs_site_counts':quotient_checks,'exact_rate_separated_generator_rows':generator_checks,
        'families':families,
        'interpretation':'Every reported completed derivative matches independent finite bit-host semantics. This is neither a general proof of the local normalizer nor a closure certificate for unfinished work.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',choices=['local','degree2','no-coop','all'],default='all')
    args=parser.parse_args();start=perf_counter();out=ROOT/'out/local-slots';out.mkdir(parents=True, exist_ok=True)
    reference=ref.check_reference()
    sources=[Path(__file__),Path(ref.__file__)]
    source_hashes={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in sources}
    reference['source_sha256']=source_hashes
    (out/'independent-reference.json').write_text(json.dumps(reference,indent=2)+'\n', encoding='utf-8')
    variants=['no-coop','local','degree2'] if args.model=='all' else [args.model]
    verified={}
    for variant in variants:
        report=verify_variant(variant);report['source_sha256']=source_hashes
        (out/variant/'verification.json').write_text(json.dumps(report,indent=2)+'\n', encoding='utf-8')
        verified[variant]={key:report[key] for key in ('passed','result_status','result_closed','independent_typed_vs_site_counts','exact_rate_separated_generator_rows')}
        print(variant,json.dumps(verified[variant]),flush=True)
    summary={'passed':True,'reference_checks':{key:value for key,value in reference.items() if key.endswith('checks') or key.endswith('counts')},
        'variants':verified,'source_sha256':source_hashes,'seconds':perf_counter()-start,
        'distinction':'Finite-host verification is separate from the explicit polynomial-degree proof of no uniform finite linear closure for positive cooperation on all cycle sizes. UNKNOWN synthesis runs remain UNKNOWN.'}
    (out/'verification-summary.json').write_text(json.dumps(summary,indent=2)+'\n', encoding='utf-8')
    print('verified in',summary['seconds'],'seconds')


if __name__=='__main__':main()
