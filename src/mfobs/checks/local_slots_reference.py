"""Independent local-slot reference; NEVER an input to composition synthesis.

Only exact finite-graph combinatorics are used. The Boolean-polynomial engine
also gives an all-cycle-size nonclosure certificate for positive cooperation.
No MolFoundry, finite-state synthesis adapter or numerical approximation is used.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction


def edges_checked(n, edges):
    assert type(n) is int and n >= 0
    result = {tuple(sorted(edge)) for edge in edges}
    assert all(len(e)==2 and 0<=e[0]<e[1]<n for e in result)
    return tuple(sorted(result))


def cycle(n):
    assert n >= 3
    return edges_checked(n, ((i,(i+1)%n) for i in range(n)))


def neighbours(n, edges):
    result=[set() for _ in range(n)]
    for u,v in edges_checked(n,edges): result[u].add(v);result[v].add(u)
    return result


def events(bits, edges, *, basal=1, cooperative=1, off=1):
    """Actual bit flips, with separate channel rates and map multiplicity."""
    n=len(bits);adj=neighbours(n,edges)
    assert all(x in (0,1) for x in bits)
    for i,x in enumerate(bits):
        changed=list(bits);changed[i]=1-x;target=tuple(changed)
        if x:
            yield target,'off',Fraction(off)
        else:
            yield target,'basal',Fraction(basal)
            # Each occupied neighbour is one distinct cooperative rule match.
            for j in sorted(adj[i]):
                if bits[j]:yield target,'cooperative',Fraction(cooperative)


def generator_value(function,bits,edges,**rates):
    before=function(bits)
    return sum(rate*(function(target)-before) for target,_,rate in events(bits,edges,**rates))


def motif_count(bits, edges, required_bits, motif_edges=(), motif_non_edges=()):
    """Injective maps of an occupancy-marked site motif into ANY finite graph.

    A None occupancy is unconstrained. Motif vertices remain distinguishable;
    automorphisms/orientations are counted, not divided out.
    """
    n=len(bits);k=len(required_bits);host=set(edges_checked(n,edges))
    positive=set(edges_checked(k,motif_edges));negative=set(edges_checked(k,motif_non_edges))
    candidates={u:tuple(v for v,x in enumerate(bits) if required_bits[u] is None or x==required_bits[u]) for u in range(k)}
    mapping={};used=set()
    def permitted(u):
        return [v for v in candidates[u] if v not in used and all(
            (tuple(sorted((u,w))) not in positive or tuple(sorted((v,z))) in host) and
            (tuple(sorted((u,w))) not in negative or tuple(sorted((v,z))) not in host)
            for w,z in mapping.items())]
    def visit():
        if len(mapping)==k:return 1
        u,choices=min(((u,permitted(u)) for u in range(k) if u not in mapping),key=lambda x:len(x[1]))
        total=0
        for v in choices:
            mapping[u]=v;used.add(v);total+=visit();used.remove(v);del mapping[u]
        return total
    return visit()


@dataclass(frozen=True)
class TypedGraph:
    nodes: tuple
    edges: tuple


def typed_slots(bits, edges):
    """Three vertices per site; None entries omit both occupancy tests."""
    nodes=tuple((3*i+j,label) for i in range(len(bits)) for j,label in enumerate(('C','P','L')))
    es=[(3*i,3*i+1,'cp') for i in range(len(bits))]+[(3*i,3*i+2,'cl') for i in range(len(bits))]
    es += [(3*i+1,3*i+2,'bound') for i,x in enumerate(bits) if x==1]
    es += [(3*u,3*v,'next') for u,v in edges_checked(len(bits),edges)]
    return TypedGraph(nodes,tuple(sorted(es)))


def typed_match_count(pattern,host,no_edges=()):
    """Independent generic typed-graph matcher; no finite slot decoder involved."""
    pl=dict(pattern.nodes);hl=dict(host.nodes)
    pe={tuple(sorted((u,v))):label for u,v,label in pattern.edges}
    he={tuple(sorted((u,v))):label for u,v,label in host.edges};ne={tuple(sorted(x)) for x in no_edges}
    incidence=lambda graph,labels:{u:Counter((edge_label,labels[b if a==u else a])
        for a,b,edge_label in graph.edges if u in (a,b)) for u in labels}
    pi=incidence(pattern,pl);hi=incidence(host,hl)
    candidates={u:tuple(v for v,label in host.nodes if label==pl[u] and
        all(hi[v][kind]>=amount for kind,amount in pi[u].items())) for u in pl}
    mapping={};used=set()
    def permitted(u):
        return [v for v in candidates[u] if v not in used and all(
            (tuple(sorted((u,w))) not in pe or he.get(tuple(sorted((v,z))))==pe[tuple(sorted((u,w)))]) and
            (tuple(sorted((u,w))) not in ne or tuple(sorted((v,z))) not in he)
            for w,z in mapping.items())]
    def visit():
        if len(mapping)==len(pl):return 1
        u,choices=min(((u,permitted(u)) for u in pl if u not in mapping),key=lambda x:len(x[1]))
        total=0
        for v in choices:
            mapping[u]=v;used.add(v);total+=visit();used.remove(v);del mapping[u]
        return total
    return visit()


def evaluate(poly, bits):
    state=sum(int(bit)<<i for i,bit in enumerate(bits))
    return sum(c for support,c in poly.items() if support & state == support)


def total_occupancy(n):
    return {1<<i:Fraction(1) for i in range(n)}


def apply_generator(poly,n,edges,*,basal=1,cooperative=1,off=1):
    """Q applied exactly to a squarefree Boolean polynomial.

    For m_S=product_{i in S}x_i, off contributes -off*|S|m_S;
    basal contributes basal*sum_i(m_{S-i}-m_S); cooperation contributes
    cooperative*sum_{i in S,j~i}(m_{(S-i)+j}-m_{S+j}).
    """
    adj=neighbours(n,edges);out=Counter()
    for support,coefficient in poly.items():
        assert 0<=support<(1<<n)
        for i in range(n):
            if not support & (1<<i):continue
            smaller=support & ~(1<<i)
            out[support]-=coefficient*Fraction(off)
            out[smaller]+=coefficient*Fraction(basal);out[support]-=coefficient*Fraction(basal)
            for j in adj[i]:
                out[smaller | (1<<j)]+=coefficient*Fraction(cooperative)
                out[support | (1<<j)]-=coefficient*Fraction(cooperative)
    return {s:c for s,c in out.items() if c}


def cycle_interval_masks(n,length):
    return {sum(1<<((start+j)%n) for j in range(length)) for start in range(n)}


PROOF = [
    'Scope: the union of all finite well-formed cycles C_n (n>=3), with Boolean occupancy and fixed nonnegative basal/off rates and fixed cooperative rate c>0.',
    'Let N=sum_i x_i. Every function on {0,1}^n has a unique multilinear polynomial representative, so differing polynomial degrees imply linear independence.',
    'For a monomial m_S, basal and off terms have degree at most |S|. Only cooperative terms -c*m_(S union {j}), with i in S and j outside S adjacent to i, can raise degree.',
    'Claim by induction: for n>=max(3,r+2), the degree-(r+1) part of Q^r N equals (-2c)^r times the sum of m_S over all proper consecutive (r+1)-site intervals S on C_n.',
    'Base r=0 is N itself. To obtain degree r+2 at the next step, extend an interval by one adjacent outside vertex. Each new proper interval has exactly two predecessor intervals, giving the next coefficient -2c times the previous coefficient. No lower-degree term can reach degree r+2 in one step.',
    'For any R, choose n>=max(3,R+2). The functions N,QN,...,Q^R N then have strictly increasing degrees 1,...,R+1 and nonzero leading coefficients. They are linearly independent on this finite Boolean cube.',
    'Consequently their restrictions cannot satisfy a common finite linear dependency on the union of all cycle sizes. There is no uniform finite-dimensional generator-invariant linear space containing N for c>0 on that language.',
    'This does not contradict fixed Ring8 closure, does not exclude approximation or specially restricted initial laws, and does not assert a generic graph-rewriting termination theorem. For c=0, QN=basal*|V|-(basal+off)*N and Q|V|=0, so a two-coordinate linear closure exists for arbitrary finite site graphs.',
]


def check_reference():
    count_checks=0;generator_checks=0;leading_checks=0
    graphs=[('C3',3,cycle(3)),('C4',4,cycle(4)),('C5',5,cycle(5)),
            ('path4',4,((0,1),(1,2),(2,3))),('star4',4,((0,1),(0,2),(0,3))),('discrete4',4,())]
    motifs=[((1,),()),((0,),()),((0,1),((0,1),)),((1,1),((0,1),)),((1,0,1),((0,1),(1,2)))]
    for _,n,edges in graphs:
        for state in range(1<<n):
            bits=tuple((state>>i)&1 for i in range(n));host=typed_slots(bits,edges)
            for required,motif_edges in motifs:
                pattern=typed_slots(required,motif_edges);ne=tuple((3*i+1,3*i+2) for i,x in enumerate(required) if x==0)
                assert motif_count(bits,edges,required,motif_edges)==typed_match_count(pattern,host,ne)
                count_checks+=1
            occupied=sum(bits);free=n-occupied;mixed=motif_count(bits,edges,(0,1),((0,1),))
            assert generator_value(sum,bits,edges)==free+mixed-occupied
            polynomial=total_occupancy(n)
            for _ in range(4):
                image=apply_generator(polynomial,n,edges)
                assert evaluate(image,bits)==generator_value(lambda x:evaluate(polynomial,x),bits,edges)
                polynomial=image;generator_checks+=1
    rows=[]
    for n in (3,4,5,6,8,12,20):
        poly=total_occupancy(n);iterations=[]
        for r in range(min(n-2,6)+1):
            degree=max(s.bit_count() for s in poly)
            leading={s:c for s,c in poly.items() if s.bit_count()==degree}
            assert degree==r+1
            assert leading=={s:Fraction((-2)**r) for s in cycle_interval_masks(n,r+1)}
            iterations.append({'iteration':r,'degree':degree,'nonzero_monomials':len(poly),
                'degree_counts':{str(k):v for k,v in sorted(Counter(s.bit_count() for s in poly).items())},
                'leading_interval_length':r+1,'leading_coefficient':str((-2)**r)})
            leading_checks+=1;poly=apply_generator(poly,n,cycle(n))
        rows.append({'cycle_size':n,'iterations':iterations})
    # c=0 is the essential contrasting case, not a hidden counterexample.
    for _,n,edges in graphs:
        image=apply_generator(total_occupancy(n),n,edges,basal=Fraction(2,5),cooperative=0,off=Fraction(3,5))
        assert image=={0:Fraction(2*n,5),**{1<<i:Fraction(-1) for i in range(n)}}
    examples=[]
    for title,bits in [('block4',(1,1,1,1,0,0,0,0)),('alternating4',(1,0,1,0,1,0,1,0))]:
        edges=cycle(8);occupied=sum(bits);free=8-occupied;mixed=motif_count(bits,edges,(0,1),((0,1),))
        examples.append({'id':title,'bits':list(bits),'occupied':occupied,'free':free,'cooperative_maps':mixed,
            'unit_rate_occupancy_derivative':free+mixed-occupied})
    assert examples[0]['occupied']==examples[1]['occupied'] and examples[0]['cooperative_maps']!=examples[1]['cooperative_maps']
    # P/L-only degree bounds: one centre can own two ports and two ligands.
    malformed=TypedGraph(((0,'C'),(1,'P'),(2,'P'),(3,'L'),(4,'L')),
                         ((0,1,'cp'),(0,2,'cp'),(0,3,'cl'),(0,4,'cl'),(1,3,'bound')))
    def pl_bounds(graph):
        labels=dict(graph.nodes)
        for u,label in graph.nodes:
            if label in ('P','L'):
                neighbour_labels=Counter(labels[b if a==u else a] for a,b,_ in graph.edges if u in (a,b))
                if neighbour_labels['C']>1 or neighbour_labels['L' if label=='P' else 'P']>1:return False
        return True
    assert pl_bounds(malformed)
    assert not pl_bounds(TypedGraph(malformed.nodes,malformed.edges+((2,3,'bound'),)))
    return {'passed':True,'independent_typed_vs_site_motif_counts':count_checks,
        'exact_generator_polynomial_checks':generator_checks,'leading_degree_checks':leading_checks,
        'occupancy_counterexample':examples,'cycle_hierarchy':rows,'all_size_proof':PROOF,
        'pl_only_bounds_counterexample':{'nodes':[list(x) for x in malformed.nodes],
            'edges':[list(x) for x in malformed.edges],'basal_adds':[2,3,'bound'],
            'NoEdge_before':True,'PL_bounds_before':True,'PL_bounds_after':False,
            'meaning':'The unchanged basal rule can bind P2 to L1, although L1 already has partner P1; C uniqueness is missing.'},
        'scope':'Independent verification only; no output here is supplied as a synthesis basis, reduction oracle or target ODE.'}


