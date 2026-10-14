"""Eight-slot semantics plug-in for the shared dynamic-basis core.

Only cycle completion and preservation evidence are scenario-specific here;
composition, derivative assembly, span extension and worklist control execute
in :mod:`mfobs.core.synthesis` and :mod:`mfobs.core.worklist`.
"""

from fractions import Fraction
import molfoundry as mf

from mfobs.core.semantic_contract import SemanticContract, machine_preservation
from mfobs.core.synthesis import Candidate, DynamicBasisAdapter, Reduction
from mfobs.core.linear_closure import RationalSpan
from mfobs.inputs import InputError
from mfobs.invariants.eight_slots import CounterFromAxioms
from mfobs.invariants.local_slots import meaning


class Adapter(DynamicBasisAdapter):
    def __init__(self, inputs):
        inputs.validate_overlap_policy()
        policy = inputs.data['composition']
        if policy['connected'] or not policy['include_empty']:
            raise InputError('All overlaps, including empty and disconnected, are required')
        if inputs.data['forbidden'] is not None or policy['mode'] not in ('public',):
            raise InputError('No external pruning theory is certified here')
        if inputs.read('reduction') != {'method': 'eight-slot-component-polynomials-v1'}:
            raise InputError('Unsupported symbolic reduction contract')
        self.inputs = inputs
        self.counter = CounterFromAxioms(inputs.read('state_language'))
        # These bounds follow from the validated Ω and hold for partial
        # patterns too. Exact degrees would incorrectly reject such patterns.
        self.composition_constraints = [
            mf.IncidenceConstraint(center, '<=', maximum,
                                   nodeLabels=neighbour, edgeLabels=edge)
            for center, neighbour, edge, maximum in (
                ('C', 'P', 'cp', 1), ('C', 'L', 'cl', 1),
                ('P', 'C', 'cp', 1), ('L', 'C', 'cl', 1),
                ('P', 'L', 'bound', 1), ('L', 'P', 'bound', 1),
                ('C', 'C', 'next', self.counter.contract['local_slots']['next_max_degree']),
            )
        ]
        self.pattern_records = {}
        self.zero = 0
        self.preservation = {}
        for name, rule in inputs.engine_rules().items():
            self.preservation[name] = self.counter.local.certify_rule_preservation(rule)
            self.preservation[name].update(
                total_C_preserved=True,
                exact_next_degree_preserved=True,
            )

        scope = (
            'all finite simple two-regular eight-C slot graphs satisfying '
            'the supplied ownership axioms'
        )
        records = tuple(
            machine_preservation(
                name,
                scope,
                certificate,
                certificate.get('passed', False)
                and certificate.get('total_C_preserved', False)
                and certificate.get('exact_next_degree_preserved', False),
            )
            for name, certificate in self.preservation.items()
        )
        contract = SemanticContract(
            domain=self.counter.contract['kind'],
            reducer='exact component-cycle completion and sparse polynomial elimination',
            exact_reduction_scope=(
                'Connected owned-slot path/cycle identity patterns completed '
                'using total_C=8 and exact next-degree two; unsupported '
                'patterns return UNKNOWN.'
            ),
            preservation=records,
        )
        super().__init__(inputs, RationalSpan(), contract)

    def reduce_pattern(self, rule, key):
        reduced = self.counter.reduce(rule)
        if reduced.status == 'zero':
            self.zero += 1

        pattern_meaning = ''
        if reduced.status == 'exact':
            pattern_meaning = meaning(reduced.local.quotient)

        return Reduction(
            status=reduced.status,
            form=reduced.polynomial,
            witness=reduced.record(),
            reason=reduced.reason,
            metadata={
                'rule': reduced.local.rule,
                'meaning': pattern_meaning,
                'record': reduced.record(),
            },
        )
    def encode_form(self, form):
        return {key: str(value) for key, value in sorted(form.items())}

    def seed_rule(self, supplied_rule, reduction):
        return reduction.metadata['rule']

    def seed_meaning(self, name, reduction):
        return reduction.metadata['meaning']

    def new_observable_name(self, state):
        return f'P{len(state.patterns)}'

    def candidate_patterns(self, rule, key, reduction):
        yield Candidate(
            key=key,
            rule=reduction.metadata['rule'],
            form=reduction.form,
            origin={'kind': 'actual-computed-jump', 'source_jump': key},
            witness=reduction.witness,
            metadata={
                'meaning': reduction.metadata['meaning'],
                'record': reduction.metadata['record'],
            },
        )

    def candidate_rank(self, candidate, residual, remainder):
        return (candidate.key,)

    def on_seed_added(self, name, reduction):
        self.pattern_records[name] = reduction.metadata['record']

    def on_candidate_added(self, name, candidate):
        self.pattern_records[name] = candidate.metadata['record']

    def finalize(self, state, composer, complete):
        return {
            'status': 'CLOSED_MODULO_EIGHT_SLOT_AXIOMS' if complete else 'UNKNOWN',
            'proof_status': 'All-overlap conditional composition with Ω-derived incidence bounds; local ownership certificates; exact cycle completion from total_C=8 and next-degree=2; sparse polynomial coefficient elimination, no state evaluation.',
            'proof_document': 'docs/eight-slot-topologies.md',
            'state_language': self.counter.contract,
            'possible_component_lengths_derived': list(self.counter.lengths),
            'rule_preservation': self.preservation,
            'symbolic_patterns': self.pattern_records, 'jump_terms': self.jump_terms,
            'basis_polynomials': {n: self.encode_form(v) for n, v in self.forms.items()},
            'symbolic_reduction_statistics': dict(self.counter.statistics),
            'zero_terms': self.zero,
            'synthesis_state_evaluations': 0, 'synthesis_orbit_evaluations': 0,
            'initial_values': {},
            'basis_kind': 'independent formal component-polynomial generators; semantic minimality not claimed',
            'limitations': [
                'The invariant theory is supplied and implemented, not automatically discovered.',
                'Disconnected site patterns and incomplete ownership remain UNKNOWN.',
                'Formal component polynomials may have additional component-inventory identities; a closed family need not be minimal.',
                'The fixed total eight is explicit global prior knowledge; no fixed topology or occupancy is supplied.',
                'The rules preserve topology: different topology sectors interact neither chemically nor probabilistically.',
            ],
        }


TOPOLOGIES = {'8': (8,), '4+4': (4, 4), '3+5': (3, 5)}


def monomial_count(key, bits, lengths):
    """Count all oriented cycle embeddings directly on the supplied host.

    Monomial 0 means no occupancy restriction, not a required free site.
    This function does not use canonicalization or completion from synthesis.
    """
    cycle, word = key.split(':'); n = int(cycle[1:])
    if len(word) != n or set(word) - {'0', '1'}:
        raise ValueError("Invalid occupancy monomial")
    ones = [i for i, bit in enumerate(word) if bit == '1']
    offset = 0; total = 0
    for size in lengths:
        if size == n:
            for start in range(n):
                for direction in (1, -1):
                    total += all(bits[offset+(start+direction*i) % n] for i in ones)
        offset += size
    return total


def polynomial_count(polynomial, bits, lengths):
    return sum(Fraction(coefficient)*monomial_count(key, bits, lengths)
               for key, coefficient in polynomial.items())
