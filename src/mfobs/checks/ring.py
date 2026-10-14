"""Re-evaluate the emitted patterns; check every rate on every finite orbit."""
import sympy as sp
from fractions import Fraction
from mfobs.core.engine import mf
from mfobs.models import ring_language as m
from mfobs.outputs import matrix

def verify(inputs,result):
    orbits,_=m.orbit_generator(8);hosts=[m.slot_ring_graph(s) for s in orbits]
    vectors=[]
    expressions = result.get('observable_expressions', {
        n: [{'coefficient': '1', 'gml': v['gml']}] for n, v in result['identity_patterns'].items()})
    evaluated = {}
    for name in result['basis']:
        vector = [Fraction()] * len(hosts)
        for term in expressions[name]:
            if term['gml'] not in evaluated:
                rule = mf.Rule.fromGMLString(term['gml'], add=False)
                evaluated[term['gml']] = m.evaluate_observable(rule, hosts)
            vector = [v + Fraction(term['coefficient']) * c for v, c in zip(vector, evaluated[term['gml']])]
        vectors.append(vector)
        if vector != list(map(Fraction, result['basis_vectors'][name])):
            raise AssertionError(f'Incorrect emitted observable {name}')
    B=sp.Matrix.hstack(*(sp.Matrix(v) for v in vectors));a=matrix(result);checks=0
    for r in result['rates']:
        rates=m.RingRates(**{q['id']:1 if q['id']==r['id'] else 0 for q in result['rates']})
        _,q=m.orbit_generator(8,rates)
        ar=a.subs({sp.Symbol(s['rate']):1 if s['id']==r['id'] else 0 for s in result['rates']})
        if q*B!=B*ar.T:raise AssertionError(f'Generator disagreement in {r["id"]}')
        checks+=B.rows*B.cols
    rank=unit_krylov_rank(result)
    if rank!=30:raise AssertionError('Unit-rate Krylov rank certificate failed')
    return {'passed':True,'unit_rate_krylov_rank':rank,'rank_certificate_prime':1000003,'orbit_count':len(orbits),'concrete_states':256,'basis_rank':B.rank(),
            'independent_pattern_recounts':B.rows*B.cols,'rate_separated_generator_checks':checks}


def unit_krylov_rank(result,prime=1000003):
    """A full rank over this finite field certifies full rank over Q."""
    a=matrix(result).subs({sp.Symbol(r['rate']):1 for r in result['rates']});n=a.rows
    vals=[[int(v.p)*pow(int(v.q),-1,prime)%prime for v in a.row(i)] for i in range(n)]
    row=[1]+[0]*(n-1);rows=[]
    for _ in range(n):
        rows.append(row);row=[sum(row[k]*vals[k][j] for k in range(n))%prime for j in range(n)]
    rank=0
    for col in range(n):
        pivot=next((i for i in range(rank,n) if rows[i][col]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inv=pow(rows[rank][col],-1,prime);rows[rank]=[(x*inv)%prime for x in rows[rank]]
        for i in range(n):
            if i!=rank:
                factor=rows[i][col];rows[i]=[(x-factor*y)%prime for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank
