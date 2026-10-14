"""Worklist adapter using ONLY local slot axioms and graph isomorphism.

No states, bit assignments, orbit counts, finite reference matrices, pattern
count vectors or precomputed word basis are imported or queried here.
"""
from collections import defaultdict
from fractions import Fraction
from mfobs.invariants.local_slots import Normalizer, meaning
from mfobs.inputs import InputError


class Adapter:
    def __init__(self, inputs):
        inputs.validate_overlap_policy()
        if inputs.data['composition']['connected'] or not inputs.data['composition']['include_empty']:
            raise InputError('This experiment requires all overlaps, including empty/disconnected ones')
        if inputs.data['forbidden'] is not None or inputs.data['composition']['mode'] not in ('public',):
            raise InputError('External overlap-pruning theories are not certified for the local-slot language')
        self.inputs = inputs
        self.normalizer = Normalizer(inputs.read('state_language'))
        config = inputs.read('reduction')
        if set(config) != {'method', 'max_pattern_sites'} or config['method'] != 'local-slot-axioms-and-isomorphism':
            raise InputError('Wrong local-slot reduction contract')
        self.max_sites = config['max_pattern_sites']
        if type(self.max_sites) is not int or self.max_sites < 1:
            raise InputError('max_pattern_sites must be positive')
        self.keys = {}; self.forms = {}; self.jump_terms = {}; self.cancelled_after_normalization = 0
        self.preservation = {name: self.normalizer.certify_rule_preservation(rule) for name, rule in inputs.engine_rules().items()}

    def initialize(self, state, composer):
        obs = self.inputs.engine_observables()
        for name in self.inputs.data['initial_observables']:
            value = self.normalizer.normalize(obs[name])
            if value.status != 'exact': raise InputError('Unsupported supplied seed: '+value.reason)
            state.add(name, value.rule, meaning(value.quotient)); self.keys[value.key] = name
            self.forms[name] = value.record()

    def expand(self, name, state, composer):
        state.equations[name] = {}
        for reaction in self.inputs.data['rules']:
            result, record = composer.pair(name, state.patterns[name].rule, reaction)
            combined = defaultdict(Fraction); reps = {}; terms = []
            for key, coefficient in sorted(result.observable_coefficients.items()):
                jump = result.observable_representatives[key]; value = self.normalizer.normalize(jump)
                self.jump_terms.setdefault(key, {'original_gml': jump.getGMLString(), 'reduction': value.record()})
                terms.append({'jump_key': key, 'coefficient': coefficient, 'status': value.status, 'reason': value.reason, 'normalized_key': value.key})
                if value.status == 'zero': continue
                if value.status == 'unknown':
                    state.unknown.append({'source': name, 'reaction': reaction['id'], 'reason': value.reason, 'jump_key': key}); continue
                combined[value.key] += Fraction(coefficient)*Fraction(reaction['per_match'])
                reps[value.key] = value
            record['terms'] = terms
            coefficients = {}
            for key, coefficient in sorted(combined.items()):
                if not coefficient:
                    self.cancelled_after_normalization += 1; continue
                value = reps[key]
                if value.quotient['site_count'] > self.max_sites:
                    state.unknown.append({'source': name, 'reaction': reaction['id'], 'reason': 'pattern site budget exceeded; not a state-language bound', 'sites': value.quotient['site_count']}); continue
                if key not in self.keys:
                    index = len(state.patterns)
                    while f'P{index}' in state.patterns:
                        index += 1
                    target = f'P{index}'; self.keys[key] = target
                    state.add(target, value.rule, meaning(value.quotient)); self.forms[target] = value.record()
                coefficients[self.keys[key]] = str(coefficient)
            state.equations[name][reaction['id']] = coefficients
            record['coefficients'] = coefficients

    def finalize(self, state, composer, complete):
        return {'status': 'CLOSED_MODULO_LOCAL_SLOT_AXIOMS' if complete else 'UNKNOWN',
                'proof_status': 'Conditional composition with all overlaps; sound local zero/equality certificates and literal condition-aware isomorphism; no finite-state basis reduction.',
                'proof_document': 'docs/local-slot-experiment.md',
                'local_contract': self.normalizer.contract, 'rule_preservation': self.preservation,
                'normalized_patterns': self.forms, 'jump_terms': self.jump_terms,
                'local_reduction_statistics': dict(self.normalizer.stats),
                'cancelled_after_local_normalization': self.cancelled_after_normalization,
                'synthesis_state_evaluations': 0, 'synthesis_orbit_evaluations': 0,
                'basis_kind': 'canonical pattern generators; independence/minimality not claimed',
                'initial_values': {},
                'limitations': ['Invariant theory is an explicitly implemented local-slot fragment, not automatically discovered or an arbitrary invariant prover.',
                    'Missing ownership in partial patterns remains UNKNOWN; it is not completed by merging vertices.',
                    'Distinct canonical patterns can still have linear counting identities; graph isomorphism is sound but not a complete linear reducer.',
                    'No fixed ring length or number of sites is supplied. A pattern-size budget is only an exploration cutoff.',
                    'Pending work or unsupported reductions mean UNKNOWN, not a closed truncated ODE.']}
