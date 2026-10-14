"""Exact injective counting on singleton components + repeated species.

Positive connected pattern components cannot cross host components. Partition
these components among the singleton hosts and copies of the repeated host.
A block placed into a single host is counted with every internal NAC and full
injective multiplicity; a NAC across distinct hosts is automatically true.
Set partitions use the least remaining bit to avoid counting blocks twice.
The result is an exact polynomial in the number m of identical host copies,
expressed in the falling-factorial basis (m)_q. NO independence assumption.
Adjacency predicates are evaluated against the COMPLETE assigned host component;
a connected host component contains all neighbours of every one of its vertices.
Hence the same partition argument applies to unary adjacency predicates. There
is no split/multiplication of correlated probabilities. This module does not
discover or validate the supplied host decomposition.
"""
from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
import molfoundry as mf
from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows, graph_key


def submasks(mask):
    """Every submask of ``mask``, from ``mask`` itself down to 0."""
    sub=mask
    while True:
        yield sub
        if sub==0:break
        sub=(sub-1)&mask


def fast_count(pattern:mf.Graph,host:mf.Graph,no_edges=(),adjacency=()):
    """Count non-induced labelled injective matches; never deduplicate maps.

    MolFoundry counts them (``Rule.numMatches``, see
    ``core.molfoundry_execution.count_matches``), always exactly.
    """
    from mfobs.core.molfoundry_execution import count_matches
    for cond in adjacency:
        if cond.center not in node_labels(pattern):raise ValueError('Adjacency center outside pattern')
    return count_matches(pattern,host,no_edges,adjacency)


def canonical_serial(pattern,no_edges,adjacency=()):
    ids={u:i for i,(u,_) in enumerate(node_rows(pattern))}
    return (make_graph({ids[u]:l for u,l in node_rows(pattern)},[(ids[u],ids[v],l) for u,v,l in edge_rows(pattern)]),tuple(sorted(tuple(sorted((ids[u],ids[v]))) for u,v in no_edges)),tuple(sorted({c.reindex(ids) for c in adjacency})))

class MixtureCounter:
    def __init__(self,repeated_host:mf.Graph,max_components=10):
        self.repeated_host=repeated_host;self.max_components=max_components
        self.count_cache={};self.polynomial_cache={}
    def single_count(self,pattern,host,no_edges,adjacency=()):
        p,nac,adj=canonical_serial(pattern,no_edges,adjacency);key=(graph_key(p),graph_key(host),nac,adj)
        if key not in self.count_cache:self.count_cache[key]=fast_count(p,host,nac,adj)
        return self.count_cache[key]
    def polynomial(self,pattern:mf.Graph,singleton_hosts:tuple[mf.Graph,...],no_edges=(),adjacency=()):
        """Return coefficients c_q such that count = sum_q c_q*(m)_q."""
        pattern,no_edges,adjacency=canonical_serial(pattern,no_edges,adjacency)
        cachekey=(graph_key(pattern),tuple(graph_key(h) for h in singleton_hosts),no_edges,adjacency)
        if cachekey in self.polynomial_cache:return self.polynomial_cache[cachekey]
        # Connected positive components (isolated vertices included).
        parent={u:u for u,_ in node_rows(pattern)}
        def root(u):
            while parent[u]!=u:u=parent[u]
            return u
        for u,v,_ in edge_rows(pattern):parent[root(u)]=root(v)
        groups=defaultdict(set)
        for u in parent:groups[root(u)].add(u)
        comps=list(groups.values());k=len(comps)
        if k>self.max_components:raise RuntimeError('too many components for exact partition counter')
        full=(1<<k)-1
        @lru_cache(None)
        def part(mask):
            ids=set().union(*(comps[j] for j in range(k) if mask>>j&1)) if mask else set()
            return (make_graph({u:l for u,l in node_rows(pattern) if u in ids},[(u,v,l) for u,v,l in edge_rows(pattern) if u in ids]), tuple((u,v) for u,v in no_edges if u in ids and v in ids),tuple(c for c in adjacency if c.center in ids))
        @lru_cache(None)
        def count_in(mask,h):
            if not mask:return 1
            p,nac,adj=part(mask)
            return self.single_count(p,h,nac,adj)
        @lru_cache(None)
        def repeat_poly(mask):
            if not mask:return (1,)
            first=mask&-mask;total=[0]*(mask.bit_count()+1)
            for sub in submasks(mask):
                if not sub&first:continue
                c=count_in(sub,self.repeated_host)
                if not c:continue
                for q,d in enumerate(repeat_poly(mask^sub)):total[q+1]+=c*d
            while len(total)>1 and not total[-1]:total.pop()
            return tuple(total)
        @lru_cache(None)
        def single_rec(mask,j):
            if j==len(singleton_hosts):return repeat_poly(mask)
            total=[0]*(k+1)
            for sub in submasks(mask):
                c=count_in(sub,singleton_hosts[j])
                if not c:continue
                for q,d in enumerate(single_rec(mask^sub,j+1)):total[q]+=c*d
            while len(total)>1 and not total[-1]:total.pop()
            return tuple(total)
        out=single_rec(full,0);self.polynomial_cache[cachekey]=out
        return out
    @staticmethod
    def evaluate(poly,m):
        v=0;ff=1
        for q,c in enumerate(poly):
            if q:ff*=m-q+1
            v+=c*ff
        return v
