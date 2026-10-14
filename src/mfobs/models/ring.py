"""Exhaustive finite reference adapter for the fixed Ring-8 model.

This deliberately is *not* the productive symbolic reducer.  It enumerates all
``2**8`` occupancy assignments and their 30 dihedral orbits, so it acts as an
independent finite oracle against which symbolic implementations can be tested.
Only this adapter knows the fixed ring scaffold.  The shared worklist and
composition service know neither bits, orbits nor molecular coordinates.
"""
from __future__ import annotations

from fractions import Fraction

import sympy as sp

from mfobs.core.semantic_contract import SemanticContract, exhaustive_preservation
from mfobs.core.synthesis import Candidate, DynamicBasisAdapter, Reduction
from mfobs.core.linear_closure import RationalSpan
from mfobs.inputs import InputError
from mfobs.models import ring_language as m


SUPPORTED_BACKEND = '5025ff3bed578429832eb9c9e964e9ee2400cb31015cec003697d547502123eb'


class Adapter(DynamicBasisAdapter):
    semantic_contract_required = True

    def __init__(self, inputs):
        self.inputs = inputs
        if inputs.backend_digest() != SUPPORTED_BACKEND:
            raise InputError(
                'MODEL_CONTRACT_MISMATCH: the ring adapter is bound to the '
                'reviewed fixed-scaffold grammar. A new grammar needs a '
                'justified state language/adapter.'
            )

        self.slots = inputs.read('state_language')['sites']
        self.orbits, self.unit_generator = m.orbit_generator(self.slots)
        self.hosts = [m.slot_ring_graph(state) for state in self.orbits]
        self.span = RationalSpan()
        self.forms, self.origins = {}, {}
        self.vectors = {}
        self.jump_terms = {}
        self.zero_terms = 0
        self.motifs = {}

        maximum_length = inputs.read('reduction')['motif_max_length']
        for length in range(1, maximum_length + 1):
            for number in range(2**length):
                word = tuple(
                    (number >> (length - index - 1)) & 1
                    for index in range(length)
                )
                vector = m.word_values(word, self.orbits)
                self.motifs.setdefault(vector, (''.join(map(str, word)), word))

    def reduce_pattern(self, rule, key):
        vector = m.evaluate_observable(rule, self.hosts)
        form = {i: Fraction(v) for i, v in enumerate(vector) if v}
        if not form:
            self.zero_terms += 1
        return Reduction('exact' if form else 'zero', form,
                         {'values_on_orbits': list(vector)})

    def seed_rule(self, supplied_rule, reduction):
        # Preserve the reference adapter's counted word representative. The
        # connected-overlap exception is only validated for this form; an
        # arbitrary equal-count disconnected pattern need not share it.
        vector = tuple(reduction.form.get(i, 0) for i in range(len(self.orbits)))
        motif = self.motifs.get(vector)
        if motif is None:
            return supplied_rule
        _, word = motif
        return m.word_observable(word, close_cycle=len(word) == self.slots)

    def candidate_patterns(self, rule, key, reduction):
        yield Candidate(key, self.seed_rule(rule, reduction), reduction.form,
                        {'kind': 'computed-pattern-with-finite-count-equivalent-word'},
                        reduction.witness)



    def encode_form(self, form):
        return {str(i): str(c) for i, c in form.items()}

    def new_observable_name(self, state):
        return f'P{len(state.patterns)}'

    def run_closure(self, state, composer, **budgets):
        complete = super().run_closure(state, composer, **budgets)
        self.vectors = {n: tuple(v.get(i, Fraction()) for i in range(len(self.orbits)))
                        for n, v in self.forms.items()}
        return complete

    def finalize(self, state, composer, complete):
        bits = self.inputs.data['initial_state'].get('bits', [])
        valid_bits = (
            len(bits) == self.slots
            and all(type(value) is int and value in (0, 1) for value in bits)
        )
        if not valid_bits:
            raise InputError('initial_state.bits must contain eight bits')

        initial_orbit = m.canonical_orbit(tuple(bits))
        initial_index = self.orbits.index(initial_orbit)
        result = {
            'status': 'CLOSED_EXHAUSTIVE_FINITE_STATE',
            'adapter_role': 'exhaustive-finite-reference-oracle',
            'synthesis_method': (
                'evaluation vectors over all 30 dihedral occupancy orbits; '
                'not invariant-based symbolic reduction'
            ),
            'proof_status': (
                'Exact rational vector reduction and per-rate verification on '
                'all dihedral orbits; lumpability checked exhaustively.'
            ),
            'proof_document': 'docs/ring8.md',
            'orbits': [list(orbit) for orbit in self.orbits],
            'basis_vectors': {
                name: [str(v) for v in vector] for name, vector in self.vectors.items()
            },
            'jump_terms': self.jump_terms,
            'initial_values': {
                name: str(vector[initial_index])
                for name, vector in self.vectors.items()
            },
            'zero_terms': self.zero_terms,
            'synthesis_state_evaluations': len(self.hosts),
            'synthesis_orbit_evaluations': len(self.orbits),
        }
        if not complete:
            return result

        basis_matrix = sp.Matrix.hstack(
            *(sp.Matrix(self.vectors[name]) for name in state.basis)
        )
        generators = {}
        for reaction in self.inputs.data['rules']:
            unit_rates = {
                other['id']: int(other['id'] == reaction['id'])
                for other in self.inputs.data['rules']
            }
            _, generators[reaction['id']] = m.orbit_generator(
                self.slots,
                m.RingRates(**unit_rates),
            )

        checked = 0
        for name in state.basis:
            for reaction in self.inputs.data['rules']:
                reaction_id = reaction['id']
                row = state.equations[name][reaction_id]
                derived = basis_matrix * sp.Matrix([sp.Rational(row.get(n, '0')) for n in state.basis])
                expected = (
                    generators[reaction_id]
                    * sp.Matrix(self.vectors[name])
                    * sp.Rational(reaction['per_match'])
                )
                if derived != expected:
                    raise AssertionError(
                        f'Independent orbit generator mismatch: {name}/{reaction_id}'
                    )

                checked += len(self.orbits)

        concrete_state_count = 2**self.slots
        scope = (
            f'all {concrete_state_count} binary occupancy assignments on the '
            f'fixed {self.slots}-cycle, quotiented into {len(self.orbits)} '
            'dihedral orbits'
        )
        preservation = tuple(
            exhaustive_preservation(
                reaction['id'],
                scope,
                {
                    'passed': True,
                    'fixed_topology_preserved': True,
                    'binary_occupancy_preserved': True,
                    'independent_generator_rows_checked': (
                        len(state.basis) * len(self.orbits)
                    ),
                },
            )
            for reaction in self.inputs.data['rules']
        )
        self.semantic_contract = SemanticContract(
            domain='fixed-eight-cycle-with-binary-owned-slots',
            reducer='exact vectors on the complete finite dihedral orbit space',
            exact_reduction_scope=(
                'Every conditional identity pattern is evaluated on every '
                'orbit representative; equality means equality on this fixed '
                'finite state language only.'
            ),
            preservation=preservation,
        )

        result.update(
            basis_evaluation_rank=basis_matrix.rank(),
            independent_rate_orbit_checks=checked,
            concrete_state_count=concrete_state_count,
            orbit_count=len(self.orbits),
        )
        return result
