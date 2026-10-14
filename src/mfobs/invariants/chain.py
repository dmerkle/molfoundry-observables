"""Exact footprint counting on a declared one-port decorated path language.

chain(0) is special. chain(n>=1) = L B^(n-1) R. Every inter-tile bond
joins the unique outgoing port of one tile to the incoming port of the next.
For a connected k-vertex pattern an embedding occupies at most k consecutive
cells. Exhausting finite *footprints* is therefore complete for all lengths.
There is no degree assumption, regression, finite-difference fit, or holdout
acceptance. Unary conditions see the complete host neighbourhood.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from dataclasses import dataclass
import hashlib,json
from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows, graph_key
from mfobs.core.conditions import matches
from mfobs.core.mixtures import canonical_serial
from mfobs.invariants.sequences import Sequence,footprint_sequence

class CountingUnknown(RuntimeError): pass

def decode(data):return make_graph(data['nodes'],data['edges'])
def components(g):
    adj={u:set() for u in node_labels(g)}
    for u,v,b in edge_rows(g):adj[u].add(v);adj[v].add(u)
    left=set(adj);out=[]
    while left:
        root=min(left);seen={root};stack=[root];left.remove(root)
        while stack:
            for v in adj[stack.pop()]&left:
                seen.add(v);stack.append(v);left.remove(v)
        out.append(seen)
    return out

def project(g,ids,no_edges=(),adjacency=()):
    return (make_graph({u:l for u,l in node_rows(g) if u in ids},[(u,v,b) for u,v,b in edge_rows(g) if u in ids and v in ids]),
            tuple((u,v) for u,v in no_edges if u in ids and v in ids),tuple(c for c in adjacency if c.center in ids))

def enumerate_maps(pattern,host,no_edges,adjacency,visit):
    """Visit every match MolFoundry admits (pattern vertex -> host vertex); return their number."""
    found=matches(pattern,host,no_edges,adjacency)
    for mapping in found:visit(mapping)
    return len(found)

class ChainLanguage:
    def __init__(self,spec):
        self.spec=spec
        if spec.get('kind')!='one-port-decorated-path' or spec.get('join_bond') not in ('-','='):
            raise ValueError('Unsupported graph-language grammar')
        self.tiles={t:decode(spec[t]) for t in ['left','repeat','right']};self.zero=decode(spec['zero'])
        self.cache={};self.bond=spec['join_bond']
        for t,g in self.tiles.items():
            for role in ['incoming','outgoing']:
                u=spec[t][role]
                if u is not None and u not in node_labels(g):raise ValueError('Port is not a tile vertex')
        if spec['left']['incoming'] is not None or spec['right']['outgoing'] is not None:
            raise ValueError('Only one-port boundaries supported')
        # Required for complete-neighbourhood unary conditions to be independent
        # of the number of repeats and the neighbouring tile's identity.
        ins=[node_labels(self.tiles[t])[spec[t]['incoming']] for t in ['repeat','right']]
        outs=[node_labels(self.tiles[t])[spec[t]['outgoing']] for t in ['left','repeat']]
        if len(set(ins))!=1 or len(set(outs))!=1:raise ValueError('Port labels must be uniform')
        self.labels=set(l for g in [self.zero,*self.tiles.values()] for u,l in node_rows(g))
        self.possible_edges=set()
        for t,g in [*self.tiles.items(),('zero',self.zero)]:
            for u,v,b in edge_rows(g):self.possible_edges.add(frozenset(((t,u),(t,v))))
        for a,b in [('left','repeat'),('repeat','repeat'),('repeat','right'),('left','right')]:
            self.possible_edges.add(frozenset(((a,spec[a]['outgoing']),(b,spec[b]['incoming']))))
        self.digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    def instantiate(self,n):
        if n<0:raise ValueError('negative length')
        if n==0:return self.zero,{u:('zero',0,u) for u,l in node_rows(self.zero)}
        if n in self.cache:return self.cache[n]
        kinds=['left']+['repeat']*(n-1)+['right'];nodes={};edges=[];locations={};prev=None
        for idx,t in enumerate(kinds):
            tile=self.tiles[t];ids={u:len(nodes)+i for i,(u,l) in enumerate(node_rows(tile))}
            nodes.update({ids[u]:l for u,l in node_rows(tile)})
            edges.extend((ids[u],ids[v],b) for u,v,b in edge_rows(tile))
            locations.update({ids[u]:(t,idx,u) for u,l in node_rows(tile)})
            inp=self.spec[t]['incoming'];out=self.spec[t]['outgoing']
            if prev is not None:edges.append((prev,ids[inp],self.bond))
            prev=ids[out] if out is not None else None
        result=make_graph(nodes,edges),locations;self.cache[n]=result;return result

@dataclass
class CountResult:
    sequence: Sequence
    witness: dict
    supports: dict

class ChainCounter:
    def __init__(self,language,joint_budget=20000):
        if type(joint_budget) is not int or joint_budget < 1:
            raise ValueError('joint_budget must be a positive integer')
        self.language=language;self.cache={};self.calls=0
        self.joint_budget=joint_budget;self.count_cache={};self._depth=0;self._joint_steps=0
    def connected(self,pattern,no_edges=(),adjacency=()):
        # Keep node IDs: supports are keyed by the actual pattern vertices.
        key=graph_key(pattern),tuple(no_edges),tuple(adjacency)
        if key in self.cache:return self.cache[key]
        self.calls+=1;k=len(node_rows(pattern))
        if not k:
            result=CountResult(Sequence.constant(1),{'method':'empty-pattern','count':1},{})
            self.cache[key]=result;return result
        if set(node_labels(pattern).values())-self.language.labels:
            result=CountResult(Sequence(),{'method':'label-absence','absent':sorted(set(node_labels(pattern).values())-self.language.labels)}, {})
            self.cache[key]=result;return result
        if len(components(pattern))!=1:raise ValueError('connected() requires a connected positive pattern')
        profiles={s:Counter() for s in ['left','right','interior','both']};supports={u:set() for u in node_labels(pattern)}
        def support(mapping,loc):
            for u,h in mapping.items():t,idx,v=loc[h];supports[u].add((t,v))
        reference_n=k+2;host,loc=self.language.instantiate(reference_n);right_idx=reference_n
        def visitor(mapping):
            indices=[loc[h][1] for h in mapping.values()];lo=min(indices);hi=max(indices)
            support(mapping,loc)
            if lo==0:
                assert hi<right_idx;profiles['left'][hi]+=1
            elif hi==right_idx:profiles['right'][right_idx-lo]+=1
            elif lo==1:profiles['interior'][hi-lo+1]+=1
        enumerate_maps(pattern,host,no_edges,adjacency,visitor)
        # Both-boundary cases occur only at n<=k-1, since all n+1 tiles
        # must be occupied by a connected embedding with k vertices.
        for n in range(1,k):
            h,ll=self.language.instantiate(n)
            def both(mapping,n=n,ll=ll):
                ids=[ll[v][1] for v in mapping.values()]
                if min(ids)==0 and max(ids)==n:
                    profiles['both'][n]+=1;support(mapping,ll)
            enumerate_maps(pattern,h,no_edges,adjacency,both)
        zero,loc0=self.language.instantiate(0)
        z=enumerate_maps(pattern,zero,no_edges,adjacency,lambda m:support(m,loc0))
        witness={'method':'connected-footprint-exhaustion','pattern_vertices':k,'reference_n':reference_n,
                 'short_both_boundary_n':list(range(1,k)),'grammar_sha256':self.language.digest,
                 **{s:dict(sorted(c.items())) for s,c in profiles.items()},'zero':z}
        result=CountResult(footprint_sequence(witness),witness,supports);self.cache[key]=result;return result
    def _spend(self):
        self._joint_steps += 1
        if self._joint_steps > self.joint_budget:
            raise CountingUnknown('joint chain counting budget exhausted; no truncated count')

    def count(self,pattern,no_edges=(),adjacency=()):
        """Exact injective counts, including disconnected and cross-NoEdge cases.

        Products count componentwise-injective maps. Their nontrivial equality
        kernels are enumerated exactly and subtracted as smaller quotient
        patterns. Cross-component NoEdge is expanded by deletion/contraction
        of the *condition*: subtract each possible positive labelled edge.
        No statistical independence or fixed polynomial degree is assumed.
        """
        no_edges=tuple(sorted(set(tuple(sorted(pair)) for pair in no_edges if pair[0]!=pair[1])))
        adjacency=tuple(sorted(set(adjacency)))
        key=graph_key(pattern),no_edges,adjacency
        if key in self.count_cache:return self.count_cache[key]
        top=self._depth==0
        if top:self._joint_steps=0
        self._depth+=1
        try:
            result=self._count(pattern,no_edges,adjacency)
            self.count_cache[key]=result
            return result
        finally:self._depth-=1

    def _count(self,pattern,no_edges,adjacency):
        cc=components(pattern)
        if len(cc)<=1:return self.connected(pattern,no_edges,adjacency)
        self._spend()
        # Unconstrained isolated vertices occupy arbitrary unused vertices of
        # their label. Every embedding of the remaining pattern uses the same
        # label inventory. This falling-factorial identity avoids Bell growth
        # for high-order moments, without a preselected moment basis.
        constrained={u for u,v,b in edge_rows(pattern)}|{v for u,v,b in edge_rows(pattern)}
        constrained|={u for pair in no_edges for u in pair}|{a.center for a in adjacency}
        free=set(node_labels(pattern))-constrained
        if free:
            core,ne,adj=project(pattern,constrained,no_edges,adjacency)
            base=self.count(core,ne,adj)
            value=base.sequence;supports=dict(base.supports);factors=[]
            used=Counter(node_labels(core).values());groups=defaultdict(list)
            for u in sorted(free):groups[node_labels(pattern)[u]].append(u)
            for label,ids in sorted(groups.items()):
                one=self.connected(make_graph({0:label}))
                for j in range(len(ids)):value=value*(one.sequence+(-used[label]-j))
                for u in ids:supports[u]=one.supports.get(0,set())
                factors.append({'label':label,'free_vertices':len(ids),'already_used':used[label],
                                'total_label_count':one.sequence.json()})
            return CountResult(value,{'method':'free-vertex-falling-factorials','core':base.witness,
                                      'factors':factors},supports)
        terms=[]
        for ids in cc:
            result=self.connected(*project(pattern,ids,no_edges,adjacency));terms.append(result)
            if not result.sequence:return CountResult(Sequence(),{'method':'impossible-connected-component','component':sorted(ids),'witness':result.witness},{})
        supports={u:ss for t in terms for u,ss in t.supports.items()}
        # A cross NoEdge is automatic for incompatible support types. Otherwise
        # its complement is the disjoint union of possible single edge labels
        # (hosts are simple). Adding that edge joins two positive components.
        owner={u:i for i,ids in enumerate(cc) for u in ids}
        for u,v in no_edges:
            if owner[u]==owner[v]:continue
            if not any(frozenset((a,b)) in self.language.possible_edges for a in supports[u] for b in supports[v]):continue
            rest=tuple(pair for pair in no_edges if pair!=(u,v))
            base=self.count(pattern,rest,adjacency);value=base.sequence;corrections=[]
            labels={b for g in [self.language.zero,*self.language.tiles.values()] for _,_,b in edge_rows(g)}
            labels.add(self.language.bond)
            for bond in sorted(labels):
                self._spend()
                augmented=make_graph(node_rows(pattern),[*edge_rows(pattern),(u,v,bond)])
                sub=self.count(augmented,rest,adjacency)
                value=value+sub.sequence*(-1)
                corrections.append({'bond':bond,'count':sub.sequence.json(),'witness':sub.witness})
            return CountResult(value,{'method':'cross-noedge-complement','endpoints':[u,v],
                                      'unguarded':base.witness,'edge_corrections':corrections},supports)
        # Retain the cheap support-disjoint product when it is justified.
        overlap=False
        for i,aa in enumerate(terms):
            for bb in terms[:i]:
                sa=set().union(*aa.supports.values());sb=set().union(*bb.supports.values())
                if sa&sb:overlap=True
        value=Sequence.constant(1)
        for t in terms:value=value*t.sequence
        if not overlap:
            return CountResult(value,{'method':'disjoint-support-product','components':[t.witness for t in terms]},supports)

        # Every componentwise-injective map has exactly one equality kernel.
        # A kernel block cannot contain two vertices of the same component.
        # Quotient injective maps biject with maps having exactly that kernel;
        # subtract each nontrivial kernel ONCE (not pairwise inclusion-exclusion).
        vertices=sorted(node_labels(pattern));blocks=[];quotients=Counter();representatives={}
        def emit():
            if len(blocks)==len(vertices):return
            mapping={u:i for i,block in enumerate(blocks) for u in block}
            try:
                graph=make_graph({i:node_labels(pattern)[block[0]] for i,block in enumerate(blocks)},
                                 [(mapping[u],mapping[v],b) for u,v,b in edge_rows(pattern)])
            except ValueError:return  # Incompatible labelled edges or loops.
            ne=tuple(sorted(set(tuple(sorted((mapping[u],mapping[v]))) for u,v in no_edges if mapping[u]!=mapping[v])))
            adj=tuple(sorted(set(a.reindex(mapping) for a in adjacency)))
            graph,ne,adj=canonical_serial(graph,ne,adj)
            key=graph_key(graph),ne,adj
            representatives.setdefault(key,graph)
            quotients[key]+=1
        def partitions(i):
            self._spend()
            if i==len(vertices):emit();return
            u=vertices[i]
            for block in blocks:
                if any(owner[v]==owner[u] or node_labels(pattern)[v]!=node_labels(pattern)[u] for v in block):continue
                common=set(supports[u])
                for v in block:common.intersection_update(supports[v])
                if not common:continue
                block.append(u);partitions(i+1);block.pop()
            blocks.append([u]);partitions(i+1);blocks.pop()
        partitions(0)
        rows=[]
        for (key,ne,adj),multiplicity in quotients.items():
            graph=representatives[key,ne,adj]
            sub=self.count(graph,ne,adj);value=value+sub.sequence*(-multiplicity)
            rows.append({'multiplicity':multiplicity,'quotient_nodes':node_rows(graph),'quotient_edges':edge_rows(graph),
                         'no_edges':ne,'adjacency':[a.gml() for a in adj],
                         'count':sub.sequence.json(),'witness':sub.witness})
        return CountResult(value,{'method':'component-kernel-partition','product_components':[t.witness for t in terms],
                                  'nontrivial_kernel_count':sum(quotients.values()),'quotient_classes':rows},supports)
