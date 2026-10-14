#!/usr/bin/env python3
"""Independent semantic checks of adjacency-condition transport and counting.

Uses small SYNTHETIC diagnostic graph rules, not additional chemical reactions.
For every fixture, directly apply the rule on complete hosts and compare the
change of a condition-bearing observable with the actual composition/jump sum.
No supplied invariant reducer, training grid, or kinetic target is used here.
"""
from mfobs.core.graph import make_graph, node_labels
from mfobs.core.rules import make_rule, rule_graph
from mfobs.core.conditions import engine_conditions, engine_no_edges
import molfoundry as mf
import time
from collections import Counter

from mfobs.core.graph import apply_match
from mfobs.core.conditions import Adjacency, matches, parse_observable, UnsupportedCondition
from mfobs.core.mixtures import fast_count,MixtureCounter
from mfobs.core.gml import disjoint_union
from mfobs.core.engine import commutator
from mfobs.models.finite_components import Counter as ComponentCounter


def G(nodes,edges=()):return make_graph(nodes,edges)

def host(d,order='-'):
    ns={0:'C',1:'O',2:'C',3:'Si'}
    es=[(0,1,'-')]
    for i in range(d):ns[4+i]='C';es.append((0,4+i,order))
    return G(ns,es)


def verify():
    start=time.perf_counter();counts=Counter();thresholds=Counter();rows=[]
    before=G({0:'C',1:'O',2:'C'},[(0,1,'-')]);after=G(node_labels(before),[(0,1,'-'),(0,2,'-')])
    dl=G(node_labels(before),[(0,1,'-'),(0,2,'=')])
    rules=[make_rule('attach',before,after,((0,2),)),
           make_rule('detach',after,before),
           make_rule('double-to-single',dl,after),
           make_rule('single-to-double',after,dl),
           make_rule('birth-center',G({3:'Si'}),G({0:'C',1:'O',2:'C',3:'Si'},[(0,1,'-'),(0,2,'-')]))]
    hosts=[host(d,order) for d in range(5) for order in ['-','=']]
    p=G({0:'C',1:'O'},[(0,1,'-')])
    for op in ['<','<=','=','>','>=']:
        for k in range(5):
            guard=Adjacency(0,('C',),('-',),op,k)
            o=make_rule('adj-observable',p,p,(),(guard,))
            erobs=o
            for r in rules:
                er=r
                c=commutator(erobs,er,connected=False,include_empty=True)
                terms=[]
                for key,coeff in c.observable_coefficients.items():
                    pp,adj=parse_observable(c.observable_representatives[key])
                    assert pp.jump_identity and not pp.generic_conditions
                    terms.append((coeff,G(pp.labels,pp.edges),pp.no_edges,adj))
                    for a in adj:thresholds[(a.op,a.count)]+=1
                for h in hosts:
                    obefore=fast_count(p,h,(),(guard,))
                    direct=0
                    for mt in matches(r, h):
                        target=apply_match(r,h,mt)
                        if target is not None:
                            direct+=fast_count(p,target,(),(guard,))-obefore
                    via=sum(co*fast_count(g,h,nac,adj) for co,g,nac,adj in terms)
                    assert direct==via,(r.name,op,k,h,direct,via)
                    counts['independent_generator_vs_composition']+=1
                rows.append({'rule':r.name,'op':op,'count':k,'hosts':len(hosts),
                   'raw_compositions':c.left_raw+c.right_raw,'jump_terms':len(terms),'passed':True})
    # First-rule condition must also survive; force different centers/thresholds.
    # Product matching is tested by the generator identity after gating r itself.
    for rk in range(1,4):
        r=make_rule(rules[0].name, rule_graph(rules[0], 'left'), rule_graph(rules[0], 'right'),
            engine_no_edges(rules[0]), (Adjacency(0,('C',),('-',),'<',rk),))
        for ok in range(1,4):
            o=make_rule('both-conditions',p,p,(),(Adjacency(0,('C',),('-',),'<',ok),))
            c=commutator(o,r)
            terms=[]
            for key,co in c.observable_coefficients.items():
                pp,adj=parse_observable(c.observable_representatives[key]);terms.append((co,G(pp.labels,pp.edges),pp.no_edges,adj))
            for h in hosts:
                val=fast_count(p,h,(),engine_conditions(o))
                direct=sum(fast_count(p,t,(),engine_conditions(o))-val
                    for mt in matches(r, h)
                    if (t:=apply_match(r,h,mt)) is not None)
                via=sum(co*fast_count(g,h,nac,adj) for co,g,nac,adj in terms)
                assert direct==via,(rk,ok,direct,via)
                counts['two_condition_transport_checks']+=1
    # Exact falling-factorial counter: constraints count all neighbours INSIDE
    # the full component to which the center is assigned, not only pattern edges.
    small=G({0:'C',1:'O',2:'C'},[(0,1,'-'),(0,2,'-')]);single=G({0:'C',1:'O'},[(0,1,'-')])
    q=disjoint_union([p,p])
    for k in [0,1,2,3]:
        adj=(Adjacency(0,('C',),('-',),'<',k),Adjacency(2,('C',),('-',),'>=',0))
        for neg in [(),((0,2),)]:
            ct=MixtureCounter(small);poly=ct.polynomial(q,(single,),neg,adj)
            for n in range(5):
                h=disjoint_union([single]+[small]*n)
                want=len(matches(q,h,neg,adj))
                assert ct.evaluate(poly,n)==fast_count(q,h,neg,adj)==want
                counts['disconnected_factorial_counter_checks']+=1
    # The stated scope is enforced: unsupported condition classes return UNKNOWN.
    bad=mf.Rule.fromGMLString('''rule [ context [node [id 0 label "C"] node [id 1 label "C"]]
       constrainShortestPath [ source 0 target 1 op "<" length 3 ] ]''',add=False)
    counter=ComponentCounter({'kind':'finite-component-multiset',
                              'components':{'C':{'nodes':[[0,'C']], 'edges':[]}}})
    try:
        counter.reduce(bad)
    except UnsupportedCondition as exc:
        res={'status':'unknown','reason':str(exc)}
    else:
        raise AssertionError('Unsupported condition accepted')
    counts['UNKNOWN_fallback_checks']+=1
    for op in ['!=','==']:
        try:Adjacency(0,('C',),('-',),op,1)
        except UnsupportedCondition:counts['rejected_unsupported_operator_checks']+=1
        else:raise AssertionError('Accepted an operator outside the engine contract')
    for nl,el in [((),('-',)),(('C',),())]:
        try:Adjacency(0,nl,el,'<',3)
        except UnsupportedCondition:counts['rejected_unsupported_filter_checks']+=1
        else:raise AssertionError('Silently accepted unsupported empty/wildcard filters')
    out={'passed':True,'counts':dict(counts),'rows':rows,
         'transported_predicate_thresholds':{f'{op} {n}':v for (op,n),v in sorted(thresholds.items())},
         'UNKNOWN_example':res,'seconds':time.perf_counter()-start,
         'scope':'finite independent diagnostic semantics, not a theorem for all nested graph conditions'}
    return out
