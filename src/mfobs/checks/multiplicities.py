"""Differential counters, including conditions and every overlap multiplicity."""
from time import perf_counter
from collections import Counter
from hashlib import sha256
import json,gc
from mfobs.core import engine as bc
from mfobs.core.forbidden import LocalPatternPredicate

def verify(inputs,result,full=True):
    rows=[];start=perf_counter();rules=inputs.engine_rules()
    expressions = result.get('observable_expressions', {
        n: [{'gml': v['gml']}] for n, v in result['identity_patterns'].items()})
    obs = {}; seen = set()
    for name, terms in expressions.items():
        for i, term in enumerate(terms):
            rule = bc.mf.Rule.fromGMLString(term['gml'], add=False)
            if rule.canonKey not in seen:
                obs[f'{name}/{i}'] = rule
                seen.add(rule.canonKey)
    spec=inputs.forbidden_spec();payload=spec.engine_payload() if spec else None
    comp=inputs.data['composition'];options={'connected':comp['connected'],'maximum':False,'include_empty':comp['include_empty']}
    def compare(name,a,b,opts):
        if payload:
            d0,c0,_=bc.grouped_restricted_common_compositions(a,b,forbidden_specification=payload,**opts)
            raw0=d0.accepted_raw
            violation=LocalPatternPredicate(payload).violation
            accepts=lambda g: violation(g) is None
            evaluator=bc.mf.RCEvaluator(intermediatePredicate=accepts,leftPredicate=accepts,rightPredicate=accepts)
        else:
            raw0,c0,_=bc.grouped_common_compositions(a,b,**opts)
            evaluator=bc.mf.RCEvaluator()
        expression=a * bc.mf.rcCommon(connected=opts['connected'],maximum=opts['maximum'],includeEmpty=opts['include_empty']) * b
        materialized=evaluator.eval(expression,onlyUnique=False)
        raw1=len(materialized);c1=Counter(r.canonKey for r in materialized)
        if c0!=c1 or raw0!=raw1:raise AssertionError(f'Multiplicity mismatch: {name}')
        rows.append({'case':name,'options':opts,'raw':raw0,'classes':len(c0),'counter_sha256':sha256(json.dumps(sorted(c0.items())).encode()).hexdigest(),'passed':True})
        print(f'  multiplicities {name}: {raw0}',flush=True);gc.collect()
    selected=list(obs) if full else list(obs)[:1]
    for o in selected:
        for rn in rules:
            for reverse in [False,True]:
                a,b=(rules[rn],obs[o]) if reverse else (obs[o],rules[rn])
                compare(f'{o}/{rn}/{"rO" if reverse else "Or"}',a,b,options)
    production=len(rows)
    return {'passed':True,'production_ordered_calls_checked':production,'comparisons':len(rows),
            'method':'production composition vs public evaluator, identical structural predicates',
            'raw_production_results':sum(r['raw'] for r in rows[:production]),'rows':rows,'seconds':perf_counter()-start}
