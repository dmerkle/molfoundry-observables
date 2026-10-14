"""Discover graph observables modulo a supplied recursive component language.

No model constructors, target equations or preselected candidate observables
are imported.
Only seed GMLs, reaction GMLs, and the declarative graph-language specification
are used. Patterns considered as representatives are projections of actual
computed jump patterns (plus the canonical empty identity).
"""
from __future__ import annotations
from mfobs.core.graph import graph_key, make_graph
from mfobs.core.rules import make_rule
from fractions import Fraction
import sympy as sp

from mfobs.core.conditions import parse_observable
from mfobs.core.synthesis import Candidate,DynamicBasisAdapter,Reduction
from mfobs.core.semantic_contract import (
    SemanticContract,
    theorem_records_from_input,
)
from mfobs.invariants.counting import Reducer,json_reduction,json_polynomial,evaluate
from mfobs.invariants.chain import components,project
from mfobs.core.linear_closure import RationalSpan as Span
from mfobs.invariants.contracts import validate_language


def expression(poly):
    n,M,x=sp.symbols('n M x',integer=True);terms=[]
    for (kind,degree,q,e),c in poly.items():
        length=n**degree if kind=='power' else sp.KroneckerDelta(n,degree)
        ff=sp.prod(M-i for i in range(q))
        terms.append(sp.Rational(c.numerator,c.denominator)*length*ff*x**e)
    return str(sp.factor(sum(terms)))

class Adapter(DynamicBasisAdapter):
    def __init__(self, inputs):
        self.inputs = inputs
        specification = inputs.read('state_language')
        configuration = inputs.read('reduction')
        inputs.validate_overlap_policy()
        if configuration.get('method') != 'exact-footprints-and-dynamic-span':
            raise ValueError('Wrong invariant discovery contract')

        self.language_contract = validate_language(
            specification,
            inputs.forbidden_spec(),
        )
        self.red = Reducer(
            specification,
            max_components=configuration.get('max_components', 10),
            joint_budget=configuration.get('joint_budget', 20_000),
        )
        self.krylov_max_dimension = configuration.get(
            'krylov_max_dimension',
            8,
        )
        valid_krylov_limit = (
            type(self.krylov_max_dimension) is int
            and 0 <= self.krylov_max_dimension <= 8
        )
        if not valid_krylov_limit:
            raise ValueError(
                'krylov_max_dimension must be an integer in 0..8; this '
                'optional symbolic diagnostic is bounded'
            )

        evidence = configuration.get('preservation_evidence')
        if evidence is None:
            raise ValueError('Recursive-language reduction needs explicit preservation_evidence')
        preservation = theorem_records_from_input(
            (rule['id'] for rule in inputs.data['rules']),
            evidence,
            {
                'definition_sha256': inputs.definition_digest(),
                'language_membership_contract': self.language_contract['status'],
                'evidence_note': (
                    'Human-readable uniform structural proof, supported by '
                    'separate finite molecular regression checks; not '
                    'machine-checked theorem proving.'
                ),
            },
        )
        contract = SemanticContract(
            domain=specification.get('kind', 'recursive-component-language'),
            reducer='exact recursive-chain footprint/partition counter',
            exact_reduction_scope=(
                'Supported conditional identity patterns over one recursive '
                'chain, one finite switch and repeated finite components; '
                'exhausted budgets return UNKNOWN.'
            ),
            preservation=preservation,
        )
        super().__init__(inputs, Span(), contract)
    def reduce_pattern(self,rule,key):
        red=self.red.reduce(key,rule)
        return Reduction(red['status'],red.get('polynomial',{}),json_reduction(red),red.get('reason',''))
    def encode_form(self,form):return json_polynomial(form)
    def meaning(self,form):return expression(form)
    def seed_origin(self,name,reduction):
        return {'kind':'input-seed','file':next(o['file'] for o in self.inputs.data['observables'] if o['id']==name),
                'counting_witness':reduction.witness}
    def _candidate(self,g,ne,adj,origin):
        key=(graph_key(g),tuple(ne),tuple(adj));red=self.red.count(g,ne,adj)
        if red['status']=='unknown' or not red['polynomial']:return None
        r=make_rule('candidate',g,g,tuple(ne),tuple(adj))
        rule=r
        return Candidate(repr(key),rule,red['polynomial'],origin,json_reduction(red),
                         {'nodes':g.numVertices,'edges':g.numEdges,'conditions':len(ne)+len(adj)})
    def after_initialize(self):
        candidate=self._candidate(make_graph({}),(),(),{'kind':'empty-identity'})
        if candidate:self.register_candidate(candidate)
    def candidate_patterns(self,rule,jump_key,reduction):
        p,adj=parse_observable(rule);g=make_graph(p.labels,p.edges);cc=components(g)
        for mask in range(1,1<<len(cc)):
            ids=set().union(*(cc[i] for i in range(len(cc)) if mask>>i&1))
            gg,ne,aa=project(g,ids,p.no_edges,adj)
            candidate=self._candidate(gg,ne,aa,{'kind':'connected-component-projection','source_jump':jump_key,'selected_vertices':sorted(ids)})
            if candidate:yield candidate
    def candidate_rank(self,candidate,residual,remainder):
        direct=set(remainder)<=set(residual);m=candidate.metadata
        return (not direct,m['nodes'],m['conditions'],m['edges'],candidate.key)
    def finalize(self,state,composer,complete):
        normal={}
        for name,p in self.forms.items():
            scale=abs(next(iter(p.values()))) if len(p)==1 else Fraction(1)
            normal[name]={'scale':str(scale),'meaning':expression({k:v/scale for k,v in p.items()})}
        initial=self.inputs.data['initial_state']
        if set(initial)!={'n','M','closed'} or any(type(initial[k]) is not int or initial[k]<0 for k in ['n','M']) or type(initial['closed']) is not bool:
            from mfobs.inputs import InputError
            raise InputError('initial_state must contain n>=0, M>=0 and closed:boolean')
        x=1-int(initial['closed'])
        out={'status':'CLOSED_MODULO_DECLARED_LANGUAGE_WITH_EXACT_COUNTING',
             'proof_status':'Exact footprint/partition counting and exact linear algebra, conditional on the supplied recursive state language and its separately justified preservation. No interpolation and no supplied finite observable basis.',
             'proof_document':'docs/chain-mixture.md','jump_terms':self.jump_terms,
             'basis_counting_functions':{k:json_polynomial(v) for k,v in self.forms.items()},
             'basis_origins':self.origins,'discovery_decisions':self.discovery_decisions,'normalisation':normal,
             'initial_values':{k:str(evaluate(p,initial['n'],initial['M'],x)) for k,p in self.forms.items()},
             'counting_statistics':{'jump_patterns':len(self.red.cache),'connected_chain_queries':self.red.chain.calls,
                                    'finite_host_counts':len(self.red.counter.finite_cache),'generated_candidates':len(self._candidates)},
             'derived_atom_caps':self.red.caps,
             'language_forbidden_contract':self.language_contract,
             'limitations':['invariants and their preservation are supplied, not automatically discovered',
                            'one-port repeated-chain grammar and two-state finite singleton supported',
                            'joint chain embeddings use exact kernel partitions; exhausted counting budgets return UNKNOWN',
                            'arbitrary nested conditions not supported; no general termination/nonclosure decision']}
        if complete and len(state.basis)<=self.krylov_max_dimension:
            from mfobs.outputs import matrix
            a=matrix({'basis':state.basis,'equations_by_reaction':state.equations,'rates':self.inputs.data['rules']})
            row=sp.zeros(1,len(state.basis));row[0,0]=1;rows=[]
            for _ in state.basis:rows.append(row);row=row*a
            out['krylov_determinant']=str(sp.factor(sp.Matrix.vstack(*rows).det()))
            out['krylov_diagnostic']={'status':'computed','scope':'optional seed-cyclicity diagnostic, not required for closure'}
        elif complete:
            out['krylov_diagnostic']={'status':'skipped','dimension':len(state.basis),
                                      'max_dimension':self.krylov_max_dimension,
                                      'reason':'optional symbolic determinant size limit; closure is fully computed, minimality not certified'}
        return out
