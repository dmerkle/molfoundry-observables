"""Exact integer-argument functions: polynomial plus finite point corrections.

Finite corrections are derived from footprint classes, never fitted to samples.
"""
from dataclasses import dataclass
from fractions import Fraction
from collections import defaultdict

@dataclass(frozen=True)
class Sequence:
    polynomial: tuple = ()
    exceptions: tuple = ()

    @classmethod
    def make(cls,polynomial=(),exceptions=()):
        p=list(map(Fraction,polynomial))
        while p and not p[-1]:p.pop()
        e={int(k):Fraction(v) for k,v in dict(exceptions).items() if v}
        return cls(tuple(p),tuple(sorted(e.items())))
    @classmethod
    def constant(cls,n):return cls.make([n])
    def __bool__(self):return bool(self.polynomial or self.exceptions)
    def value(self,n):
        if n<0:raise ValueError('n must be nonnegative')
        return sum(c*n**i for i,c in enumerate(self.polynomial))+dict(self.exceptions).get(n,0)
    def __add__(self,b):
        if not isinstance(b,Sequence):b=Sequence.constant(b)
        p=[Fraction(0)]*max(len(self.polynomial),len(b.polynomial));e=defaultdict(Fraction,dict(self.exceptions))
        for i,c in enumerate(self.polynomial):p[i]+=c
        for i,c in enumerate(b.polynomial):p[i]+=c
        for k,c in b.exceptions:e[k]+=c
        return Sequence.make(p,e)
    __radd__=__add__
    def __mul__(self,b):
        if not isinstance(b,Sequence):b=Sequence.constant(b)
        p=[Fraction(0)]*max(0,len(self.polynomial)+len(b.polynomial)-1)
        for i,c in enumerate(self.polynomial):
            for j,d in enumerate(b.polynomial):p[i+j]+=c*d
        e={}
        for n in dict(self.exceptions).keys()|dict(b.exceptions).keys():
            e[n]=self.value(n)*b.value(n)-sum(c*n**i for i,c in enumerate(p))
        return Sequence.make(p,e)
    __rmul__=__mul__
    def json(self):return {'polynomial':[str(x) for x in self.polynomial], 'point_corrections':{str(k):str(v) for k,v in self.exceptions}}


def footprint_sequence(witness):
    """Replay a connected-embedding footprint count for every n>=0."""
    p=[Fraction(0),Fraction(0)];e=defaultdict(Fraction)
    for side in ['left','right']:
        for j,count in witness[side].items():
            j=int(j);count=int(count);p[0]+=count
            for n in range(j+1):e[n]-=count
    for width,count in witness['interior'].items():
        w=int(width);count=int(count);p[0]-=count*w;p[1]+=count
        for n in range(w+1):e[n]-=count*(n-w)
    for n,count in witness['both'].items():e[int(n)]+=int(count)
    e[0]+=int(witness['zero'])
    return Sequence.make(p,e)
