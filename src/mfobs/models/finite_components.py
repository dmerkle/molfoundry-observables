"""Exact connected-observable synthesis on arbitrary finite-species multisets.

This deliberately narrow baseline counts no disconnected nonzero patterns.
The supplied component library is not a target observable basis; representatives
and span extensions come only from actual computed commutator jump inputs.
"""
from mfobs.core.graph import edge_rows, graph_key, has_edge, make_graph, node_labels, node_rows
from mfobs.core.rules import preserves_vertices
from fractions import Fraction
import re

from mfobs.core.conditions import UnsupportedCondition, identity_rule, matches, parse_observable
from mfobs.core.graph import apply_match
from mfobs.core.mixtures import canonical_serial, fast_count
from mfobs.inputs import InputError
from mfobs.invariants.chain import components
from mfobs.invariants.contracts import _graph
from mfobs.core.linear_closure import RationalSpan as Span

ONE = '__ONE__'


def encoded(poly):
    return {key: str(value) for key, value in sorted(poly.items()) if value}


def expression(poly):
    return ' + '.join(f'({value})*{key}' for key, value in sorted(poly.items())) or '0'


def identity(graph, no_edges=(), adjacency=()):
    return identity_rule(graph, 'count', tuple(no_edges), tuple(adjacency))


class Counter:
    def __init__(self, language):
        if (set(language) != {'kind', 'components'} or
                language['kind'] != 'finite-component-multiset' or
                not isinstance(language['components'], dict) or not language['components']):
            raise InputError('Expected nonempty finite-component-multiset language')
        self.hosts = {}; seen = set(); self.cache = {}
        for name, data in sorted(language['components'].items()):
            if (not isinstance(name, str) or name == ONE or
                    not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', name)):
                raise InputError('Component names must be identifiers, excluding '+ONE)
            graph = _graph(data, name)
            if len(components(graph)) != 1:
                raise InputError('Every library component must be nonempty and connected')
            key = identity(graph).canonKey
            if key in seen:
                raise InputError('Library components must be pairwise nonisomorphic')
            seen.add(key); self.hosts[name] = graph

    def count(self, graph, no_edges=(), adjacency=()):
        canonical, ne, adj = canonical_serial(graph, no_edges, adjacency)
        key = (graph_key(canonical), ne, adj)
        if key not in self.cache:
            try:
                g = canonical
                if any(has_edge(g, u, v) for u, v in ne):
                    poly = {}
                elif not node_rows(g):
                    poly = {ONE: Fraction(1)}
                elif len(components(g)) != 1:
                    raise UnsupportedCondition('Disconnected noncancelled jump: multivariate counting not implemented')
                else:
                    poly = {name: Fraction(fast_count(g, host, ne, adj))
                            for name, host in self.hosts.items()}
                    poly = {name: value for name, value in poly.items() if value}
                result = {'status': 'exact' if poly else 'zero', 'polynomial': poly,
                          'method': 'connected-input-count-times-independent-copy-number'}
            except UnsupportedCondition as exc:
                result = {'status': 'unknown', 'reason': str(exc), 'polynomial': {}}
            self.cache[key] = result
        return self.cache[key]

    def reduce(self, rule):
        parsed, adjacency = parse_observable(rule)
        if (not parsed.jump_identity or parsed.generic_conditions or
                not preserves_vertices(rule)):
            raise UnsupportedCondition('Only vertex-preserving conditional identity jumps are supported')
        return self.count(make_graph(parsed.labels, parsed.edges), parsed.no_edges, adjacency)

    def preservation(self, rules):
        """Connected L localizes every application to exactly one component."""
        records = []
        for name, rule in rules.items():
            parsed, adjacency = parse_observable(rule)
            if not preserves_vertices(rule):
                raise InputError('Reactions must preserve every vertex and its label')
            left = make_graph(parsed.labels, parsed.edges)
            if len(components(left)) != 1:
                raise InputError('Reaction inputs must be nonempty and connected')
            oracle = rule
            for species, host in self.hosts.items():
                maps = matches(left, host, parsed.no_edges, adjacency)
                valid = 0
                for mapping in maps:
                    target = apply_match(oracle, host, mapping)
                    if target is None:
                        continue
                    for ids in components(target):
                        part = make_graph({u: node_labels(target)[u] for u in ids},
                            [(u, v, bond) for u, v, bond in edge_rows(target) if u in ids])
                        if not any(part.numVertices == h.numVertices and part.numEdges == h.numEdges
                            and fast_count(part, h)
                            for h in self.hosts.values()):
                            raise InputError(f'Language not preserved: {name} on {species}')
                    valid += 1
                records.append({'reaction': name, 'component': species,
                                'raw_input_matches': len(maps), 'admissible_applications': valid})
        return {'status': 'VERIFIED_ALL_MULTISETS', 'checks': records,
                'argument': 'Connected nonempty L lies in one host component; exhaustive local successors split only into declared components. Arbitrary remaining components are unchanged.'}


from mfobs.core.synthesis import DynamicBasisAdapter, Reduction


class Adapter(DynamicBasisAdapter):
    semantic_contract_required = False
    def __init__(self, inputs):
        self.inputs = inputs; config = inputs.read('reduction')
        if (set(config)-{'method'} or
                config.get('method') != 'exact-connected-components-and-dynamic-span'):
            raise InputError('Wrong connected finite-component reduction contract')
        if inputs.data['composition']['mode'].startswith('restricted') or inputs.forbidden_spec():
            raise InputError('This baseline requires unrestricted composition and no pruning theory')
        self.counter = Counter(inputs.read('state_language'))
        self.contract = self.counter.preservation(inputs.engine_rules())
        self.span = Span(); self.forms = {}; self.origins = {}; self.jump_terms = {}

    def reduce_pattern(self, rule, key):
        try:
            red = self.counter.reduce(rule)
        except (UnsupportedCondition, ValueError) as exc:
            return Reduction('unknown', {}, reason=str(exc))
        witness = {**red, 'polynomial': encoded(red['polynomial'])}
        return Reduction(red['status'], red['polynomial'], witness, red.get('reason', ''))


    def encode_form(self, form):
        return encoded(form)

    def meaning(self, form):
        return expression(form)

    def finalize(self, state, composer, complete):
        initial = self.inputs.data['initial_state']
        populations = initial.get('populations', {})
        if (set(initial) != {'populations'} or set(populations) != set(self.counter.hosts) or
                any(type(v) is not int or v < 0 for v in populations.values())):
            raise InputError('initial_state requires one nonnegative integer population per species')
        values = {**populations, ONE: 1}
        return {'status': 'CLOSED_MODULO_FINITE_COMPONENT_MULTISET_LANGUAGE',
                'proof_status': 'Exact connected-pattern counting and rational span closure; connected-rule preservation verified on the complete finite component library, for arbitrary multiplicities.',
                'basis_counting_functions': {k: encoded(v) for k, v in self.forms.items()},
                'basis_origins': self.origins, 'jump_terms': self.jump_terms,
                'preservation_certificate': self.contract,
                'initial_values': {k: str(sum(c*values[s] for s, c in p.items())) for k, p in self.forms.items()},
                'normalisation': {k: {'scale': '1', 'meaning': expression(p)} for k, p in self.forms.items()},
                'limitations': ['Finite component library is supplied, not discovered; arbitrary independent copy numbers.',
                                'Disconnected noncancelled jump patterns or exhausted counting budgets give UNKNOWN.',
                                'Only NoEdge and supported unary adjacency conditions; no general NAC separation claim.']}
