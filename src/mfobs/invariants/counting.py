"""Exact counts on a declared chain + one finite switch + repeated molecules.

The coordinate algebra has n>=0, M>=0 and x in {0,1}. It has no supplied
finite target basis and no bound on monomer degree. Point corrections in n
are retained rather than fitted away. Budgets and unsupported cases fail closed.
"""
from __future__ import annotations
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from mfobs.core.graph import make_graph, node_labels, graph_key
from mfobs.core.conditions import parse_observable,UnsupportedCondition
from mfobs.core.mixtures import canonical_serial,fast_count,submasks
from mfobs.invariants.chain import ChainLanguage,ChainCounter,decode,components,project,CountingUnknown
from mfobs.invariants.sequences import Sequence

class ExactMixtureCounter:
    def __init__(self,chain,switch_states,repeated,*,max_components=10):
        self.chain=chain;self.switch_states=switch_states;self.repeated=repeated;self.max_components=max_components
        self.finite_cache={};self.cache={};self.witnesses={}
    def finite(self,p,h,ne=(),adj=()):
        p,ne,adj=canonical_serial(p,ne,adj);key=graph_key(p),graph_key(h),ne,adj
        if key not in self.finite_cache:self.finite_cache[key]=fast_count(p,h,ne,adj)
        return self.finite_cache[key]
    def polynomial(self,p,no_edges=(),adjacency=()):
        p,no_edges,adjacency=canonical_serial(p,no_edges,adjacency);key=graph_key(p),no_edges,adjacency
        if key in self.cache:return self.cache[key]
        cc=components(p);k=len(cc)
        if k>self.max_components:raise CountingUnknown('component partition budget exceeded')
        full=(1<<k)-1
        @lru_cache(None)
        def part(mask):
            ids=set().union(*(cc[i] for i in range(k) if mask>>i&1))
            return project(p,ids,no_edges,adjacency)
        compatible={'chain':0,'monomer':0,'open':0,'closed':0}
        component_witness=[]
        for i in range(k):
            g,ne,adj=part(1<<i);ch=self.chain.connected(g,ne,adj)
            row={'ids':sorted(cc[i]),'chain':ch.sequence.json(),'chain_witness':ch.witness}
            if ch.sequence:compatible['chain']|=1<<i
            for name,h in [('monomer',self.repeated),*self.switch_states.items()]:
                c=self.finite(g,h,ne,adj);row[name]=c
                if c:compatible[name]|=1<<i
            component_witness.append(row)
            if not any(mask>>i&1 for mask in compatible.values()):
                result={};proof={'method':'no-host-for-connected-component','component':row}
                self.cache[key]=result;self.witnesses[key]=proof;return result
        @lru_cache(None)
        def repeat(mask):
            if mask & ~compatible['monomer']:return ()
            if not mask:return (1,)
            first=mask&-mask;total=[0]*(mask.bit_count()+1)
            for sub in submasks(mask):
                if not sub&first:continue
                c=self.finite(*[part(sub)[0],self.repeated,*part(sub)[1:]])
                if not c:continue
                for q,d in enumerate(repeat(mask^sub)):total[q+1]+=c*d
            while total and not total[-1]:total.pop()
            return tuple(total)
        @lru_cache(None)
        def chain_count(mask):
            if mask & ~compatible['chain']:return Sequence()
            return self.chain.count(*part(mask)).sequence
        @lru_cache(None)
        def switch_count(mask,name):
            if mask & ~compatible[name]:return 0
            g,ne,adj=part(mask);return self.finite(g,self.switch_states[name],ne,adj)
        all_values={};allocations={}
        for switch in ['closed','open']:
            vals=[Sequence()]*(k+1);proof=[]
            for catmask in submasks(full & compatible[switch]):
                catcount=switch_count(catmask,switch)
                if not catcount:continue
                remain=full^catmask
                for chmask in submasks(remain & compatible['chain']):
                    rest=remain^chmask;rp=repeat(rest)
                    if not rp:continue
                    ch=chain_count(chmask)
                    if not ch:continue
                    for q,c in enumerate(rp):
                        if c:vals[q]=vals[q]+ch*(catcount*c)
                    proof.append({'switch_mask':catmask,'chain_mask':chmask,'repeat_mask':rest,
                                  'switch_count':catcount,'chain_sequence':ch.json(),'repeat_factorial_coefficients':list(rp)})
            all_values[switch]=vals;allocations[switch]=proof
        poly={}
        def emit(seq,q,e):
            for degree,c in enumerate(seq.polynomial):
                if c:poly['power',degree,q,e]=c
            for n,c in seq.exceptions:
                if c:poly['point',n,q,e]=c
        for q,(cl,op) in enumerate(zip(all_values['closed'],all_values['open'])):
            emit(cl,q,0);emit(op+cl*(-1),q,1)
        self.cache[key]=poly
        self.witnesses[key]={'method':'component-partition-falling-factorial','components':component_witness,'allocations':allocations,
                             'boolean_quotient':'f = f_closed + x*(f_open-f_closed); x^2=x'}
        return poly

class Reducer:
    def __init__(self,language,*,max_components=10,joint_budget=20000):
        self.language=language;self.chain=ChainCounter(ChainLanguage(language['chain']),joint_budget=joint_budget)
        self.counter=ExactMixtureCounter(self.chain,{s:decode(language[s]) for s in ['open','closed']},decode(language['monomer']),max_components=max_components)
        self.cache={};self.pattern_cache={}
        # Derive finite atom caps from the supplied component grammar, not from
        # observable names or stored answers. Repeated species make caps infinite.
        b=Counter(node_labels(self.chain.language.tiles['repeat']).values());rep=Counter(node_labels(self.counter.repeated).values())
        labels=set(self.chain.language.labels)|set(node_labels(self.counter.repeated).values())
        for g in self.counter.switch_states.values():labels.update(node_labels(g).values())
        self.caps={}
        for l in labels:
            if b[l] or rep[l]:continue
            cap=max(Counter(node_labels(self.chain.language.zero).values())[l],sum(Counter(node_labels(self.chain.language.tiles[t]).values())[l] for t in ['left','right']))
            cap+=max(Counter(node_labels(g).values())[l] for g in self.counter.switch_states.values())
            self.caps[l]=cap
    def count(self,p,no_edges=(),adjacency=()):
        p,no_edges,adjacency=canonical_serial(p,no_edges,adjacency);key=graph_key(p),no_edges,adjacency
        if key in self.pattern_cache:return self.pattern_cache[key]
        for l,n in Counter(node_labels(p).values()).items():
            if l in self.caps and n>self.caps[l]:
                result={'status':'zero','polynomial':{},'witness':{'method':'derived-atom-cap','label':l,'cap':self.caps[l],'requested':n}}
                self.pattern_cache[key]=result;return result
        try:
            poly=self.counter.polynomial(p,no_edges,adjacency)
            result={'status':'exact' if poly else 'zero','polynomial':poly,'witness':self.counter.witnesses[key]}
        except CountingUnknown as exc:
            result={'status':'unknown','reason':str(exc),'polynomial':{}}
        self.pattern_cache[key]=result;return result
    def reduce(self,key,rule):
        if key in self.cache:return self.cache[key]
        try:
            p,adj=parse_observable(rule)
            if p.generic_conditions or not p.jump_identity:raise UnsupportedCondition('unsupported conditional jump span')
            from mfobs.core.rules import preserves_vertices
            if not preserves_vertices(rule):
                raise UnsupportedCondition('jump with vertex deletion/recreation needs a dangling-aware counter')
            result=self.count(make_graph(p.labels,p.edges),p.no_edges,adj)
        except (UnsupportedCondition,ValueError) as exc:result={'status':'unknown','reason':str(exc),'polynomial':{}}
        self.cache[key]=result;return result

def json_polynomial(poly):return {'|'.join(map(str,k)):str(v) for k,v in sorted(poly.items()) if v}
def json_reduction(red):return {**{k:v for k,v in red.items() if k!='polynomial'},'polynomial':json_polynomial(red.get('polynomial',{}))}

def evaluate(poly,n,m,x):
    out=Fraction()
    for (kind,degree,q,e),c in poly.items():
        ff=1
        for i in range(q):ff*=m-i
        out+=c*(n**degree if kind=='power' else int(n==degree))*ff*x**e
    return out
