"""Symbolic component-polynomial reduction from size/degree axioms.

No host graphs, occupancy assignments, orbit tables or evaluation matrices are
used. Monomials are generated on demand from actual composed count patterns.
Binary words here mark polynomial factors x_i, NOT full occupancy states:
0 means unspecified, while a free-site test is expanded as (1-x_i).
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

from mfobs.invariants.local_slots import CONTRACT as LOCAL_CONTRACT, Normalizer


CONTRACT = {
    'kind': 'eight-two-regular-owned-slots-v1',
    'local_slots': {**LOCAL_CONTRACT, 'next_max_degree': 2},
    'total_C': 8,
    'each_C_exactly_two_next_neighbours': True,
    'require_connected': False,
}


def canonical_monomial(length, ones):
    word = ''.join('1' if i in ones else '0' for i in range(length))
    reverse = word[::-1]
    canonical = min([word[i:]+word[:i] for i in range(length)] +
                    [reverse[i:]+reverse[:i] for i in range(length)])
    return f'C{length}:{canonical}'


def word_polynomial(length, word):
    """Raw cycle-embedding count with free tests expanded by inclusion-exclusion."""
    fixed = {i for i, state in enumerate(word) if state == 1}
    free = [i for i, state in enumerate(word) if state == 0]
    result = Counter()
    for size in range(len(free)+1):
        for subset in combinations(free, size):
            result[canonical_monomial(length, fixed | set(subset))] += (-1)**size
    return {key: Fraction(value) for key, value in result.items() if value}


@dataclass
class PolynomialReduction:
    status: str
    reason: str
    local: object
    polynomial: dict
    completions: tuple = ()

    def record(self):
        return {'status': self.status, 'reason': self.reason,
                'polynomial': {key: str(value) for key, value in sorted(self.polynomial.items())},
                'completion_certificates': list(self.completions),
                'local_reduction': self.local.record()}


class CounterFromAxioms:
    def __init__(self, contract):
        if contract != CONTRACT:
            raise ValueError('Unsupported eight-slot invariant contract')
        self.contract = contract
        self.local = Normalizer(contract['local_slots'])
        self.size = contract['total_C']
        # A finite simple 2-regular graph is a union of cycles >=3. A component
        # of size m fits iff the remainder is zero or can itself be a cycle.
        self.lengths = tuple(m for m in range(3, self.size+1)
                             if self.size-m == 0 or self.size-m >= 3)
        self.cache = {}
        self.statistics = Counter()

    def reduce(self, rule):
        key = rule.canonKey
        if key not in self.cache:
            local = self.local.normalize(rule)
            result = self._reduce_local(local)
            self.cache[key] = result
            self.statistics[result.status] += 1
        return self.cache[key]

    def _reduce_local(self, local):
        def done(status, reason, polynomial=None, completions=()):
            return PolynomialReduction(status, reason, local, polynomial or {}, tuple(completions))
        if local.status != 'exact':
            return done(local.status, local.reason)
        q = local.quotient
        n = q['site_count']
        if n > self.size:
            return done('zero', 'More pattern C vertices than the declared total_C=8; injectivity forbids a match')
        if n == 0:
            # H_m(empty)=2m * number of m-cycle components. Sum_m H_m=2*8.
            return done('exact', 'Empty identity via the supplied total number of C vertices',
                        {canonical_monomial(m, set()): Fraction(1, 2*self.size) for m in self.lengths},
                        [{'argument': 'sum_m H_m(empty)=2*total_C'}])
        adjacent = {i: set() for i in range(n)}
        for a, b in q['next_edges']:
            adjacent[a].add(b); adjacent[b].add(a)
        reached = {0}; pending = [0]
        while pending:
            for v in adjacent[pending.pop()]-reached:
                reached.add(v); pending.append(v)
        if len(reached) != n:
            return done('unknown', 'Disconnected site patterns are outside this component-polynomial fragment')
        ends = [i for i in adjacent if len(adjacent[i]) == 1]
        is_cycle = n >= 3 and all(len(v) == 2 for v in adjacent.values())
        if not is_cycle and n != 1 and len(ends) != 2:
            return done('unknown', 'Pattern is neither a simple path nor a cycle')
        start = min(adjacent) if is_cycle or n == 1 else min(ends)
        order = [start]
        while len(order) < n:
            following = sorted(adjacent[order[-1]]-set(order))
            if not following:
                return done('unknown', 'Could not order the connected component')
            order.append(following[0])
        positions = {old: i for i, old in enumerate(order)}
        word = [q['occupancies'][i] for i in order]
        lengths = (n,) if is_cycle and n in self.lengths else () if is_cycle else self.lengths
        total = Counter(); certificates = []
        for m in lengths:
            if m < n:
                continue
            cycle_edges = {tuple(sorted((i, (i+1) % m))) for i in range(m)}
            if any(tuple(sorted((positions[a], positions[b]))) in cycle_edges
                   for a, b in q['no_next_edges']):
                continue
            # A path with a specified first edge has a unique extension around
            # its host cycle. At a single C there are two choices of direction.
            factor = Fraction(1, 2) if n == 1 else Fraction(1)
            for key, coefficient in word_polynomial(m, word).items():
                total[key] += factor*coefficient
            certificates.append({'cycle_length': m, 'pattern_sites': n,
                'extension_multiplicity': 2 if n == 1 else 1,
                'count_factor': str(factor), 'ordered_original_sites': order,
                'occupancy_tests': word,
                'added_sites_are_unspecified': m-n,
                'argument': 'Exact degree two completes a connected match to its unique cycle component; total_C rules out a remainder of one or two vertices.'})
        polynomial = {key: coefficient for key, coefficient in total.items() if coefficient}
        if not polynomial:
            return done('zero', 'No cycle completion satisfies size, degree and NoEdge constraints', completions=certificates)
        return done('exact', 'Cycle-completion bijection and polynomial expansion of free-site NoEdges', polynomial, certificates)
