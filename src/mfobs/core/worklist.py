"""Shared FIFO closure driver and rate-separated composition service.

Adapters supply state-language reduction. No chemical model or finite-state
enumeration is built into this module. Hitting an exploration limit is UNKNOWN,
never a non-closure proof or a truncated closed ODE.
"""
from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from time import perf_counter

from mfobs.core import engine as bc

@dataclass
class Pattern:
    name: str
    rule: object
    depth: int
    meaning: str = ''

class State:
    def __init__(self, adapter, progress=False, checkpoint=None):
        self.patterns = {}
        self.queue = deque()
        self.steps = []
        self.equations = {}
        self.unknown = []
        self.adapter = adapter
        self.progress = progress
        self.checkpoint = checkpoint
        self.current_depth = -1

    def record_expansion(self, step):
        self.steps.append(step)
        if self.progress:
            print(f"expanded {step['expanded']}: basis={len(self.patterns)}, "
                  f"queue={len(step['queue'])}, unknown={len(self.unknown)}", flush=True)
        if self.checkpoint:
            partial = {'closed': False, 'partial': True, 'basis': self.basis,
                       'steps': self.steps, 'pending': step['queue'], 'unknown': self.unknown}
            Path(self.checkpoint).write_text(json.dumps(partial, indent=2) + '\n', encoding='utf-8')

    def add(self, name, rule, meaning=''):
        if name not in self.patterns:
            self.patterns[name] = Pattern(
                name,
                rule,
                self.current_depth + 1,
                meaning,
            )
            self.queue.append(name)
        return name

    @property
    def basis(self):
        return list(self.patterns)

class Composer:
    def __init__(self, inputs, progress=False, record_compositions=False, constraints=None):
        self.constraints = constraints
        self.record_compositions = record_compositions
        self.inputs = inputs
        self.rules = inputs.engine_rules()
        self.records = []
        self.progress = progress
        self.settings = inputs.data['composition']
        self.spec = inputs.forbidden_spec()
        self.payload = (
            self.spec.engine_payload()
            if self.settings['mode'].startswith('restricted')
            else None
        )

    def pair(self, source, observable, reaction):
        reaction_id = reaction['id']
        if self.progress:
            print(f'  compose {source} / {reaction_id}', flush=True)

        start = perf_counter()
        commutator = bc.commutator(
            observable,
            self.rules[reaction_id],
            connected=self.settings['connected'],
            include_empty=self.settings['include_empty'],
            forbidden_specification=self.payload,
            constraints=self.constraints,
        )
        record = {
            'source': source,
            'reaction': reaction_id,
            'left_raw': commutator.left_raw,
            'right_raw': commutator.right_raw,
            'seconds': perf_counter() - start,
            'left_restricted': (
                asdict(commutator.left_restricted)
                if commutator.left_restricted
                else None
            ),
            'right_restricted': (
                asdict(commutator.right_restricted)
                if commutator.right_restricted
                else None
            ),
            'rate_per_match_multiplier': reaction['per_match'],
        }
        if self.record_compositions:
            # Serialize before a backend uses a locality certificate to reduce
            # the sum. Include canceled conditional classes with coefficient 0.
            record['conditional_composites'] = [
                {'key': key, 'coefficient': commutator.composite_coefficients.get(key, 0),
                 'gml': rule.getGMLString()}
                for key, rule in sorted(commutator.composite_representatives.items())
            ]
            record['projected_terms'] = [
                {'key': key, 'coefficient': commutator.observable_coefficients.get(key, 0),
                 'gml': rule.getGMLString()}
                for key, rule in sorted(commutator.observable_representatives.items())
            ]
        self.records.append(record)
        return commutator, record


def _semantic_contract(adapter, state):
    """Return the adapter evidence and fail closed when it is insufficient."""
    if not hasattr(adapter, 'semantic_contract'):
        if getattr(adapter, 'semantic_contract_required', False):
            state.unknown.append({
                'reason': 'required semantic contract was not produced'
            })
        return None
    contract = adapter.semantic_contract.record()
    if not contract['allows_conditional_closure']:
        state.unknown.append({
            'reason': (
                'rule preservation contract does not allow conditional closure'
            )
        })
    return contract


def run(inputs, *, progress=False, max_expansions=1000, checkpoint=None, record_compositions=False, maximum_dimension=1000):
    """Execute composition/reduction until closure or an honest UNKNOWN.

    ``record_compositions`` serializes conditional classes and projected terms
    before backend reduction for inspection; it does not change the worklist.
    """
    if max_expansions < 0:
        raise ValueError('max_expansions must be nonnegative')

    mode = inputs.data['composition']['mode']

    start = perf_counter()
    adapter = inputs.adapter()
    composer = Composer(inputs, progress, record_compositions=record_compositions,
                        constraints=getattr(adapter, 'composition_constraints', None))
    state = State(adapter, progress, checkpoint)
    adapter.initialize(state, composer)
    limited = False

    if hasattr(adapter, 'run_closure'):
        limited = not adapter.run_closure(state, composer, max_expansions=max_expansions,
                                           maximum_dimension=maximum_dimension)

    while state.queue and not hasattr(adapter, 'run_closure'):
        if len(state.steps) >= max_expansions:
            limited = True
            state.unknown.append({
                'reason': 'expansion budget exhausted',
                'pending': list(state.queue),
            })
            break

        name = state.queue.popleft()
        state.current_depth = state.patterns[name].depth
        before = state.basis
        adapter.expand(name, state, composer)
        state.record_expansion({
            'expanded': name,
            'before': before,
            'after': state.basis,
            'queue': list(state.queue),
        })

    algorithm_complete = not limited and not state.unknown
    metadata = adapter.finalize(
        state,
        composer,
        complete=algorithm_complete,
    )
    semantic_contract = _semantic_contract(adapter, state)
    closed = not limited and not state.unknown

    result = {
        'schema_version': 1,
        'model': inputs.data['id'],
        'title': inputs.data['title'],
        'closed': closed,
        'basis': state.basis,
        'initial_observables': list(inputs.data['initial_observables']),
        'depths': {
            name: pattern.depth for name, pattern in state.patterns.items()
        },
        'steps': state.steps,
        'equations_by_reaction': state.equations,
        'unknown': state.unknown,
        'pending': list(state.queue),
        'raw_compositions': sum(
            record['left_raw'] + record['right_raw']
            for record in composer.records
        ),
        'ordered_composition_calls': 2 * len(composer.records),
        'composition_mode': mode,
        'expansions': composer.records,
        'normalisation': inputs.data['normalisation'],
        'input_initial_state': inputs.data['initial_state'],
        'rates': [
            {
                key: reaction[key]
                for key in ['id', 'rate', 'per_match', 'kind']
            }
            for reaction in inputs.data['rules']
        ],
        'input_definition_sha256': inputs.definition_digest(),
        'identity_patterns': {
            name: {
                'gml': pattern.rule.getGMLString(),
                'meaning': pattern.meaning,
            }
            for name, pattern in state.patterns.items() if pattern.rule is not None
        },
        'forbidden_theory': (
            composer.spec.certificate() if composer.payload else None
        ),
        'specification_layers': {
            'chemistry_and_domain': inputs.specification_certificate(),
            'initial_mixture': inputs.initial_mixture(),
            'initial_mixture_used_for_symbolic_closure': False,
        },
        'semantic_contract': semantic_contract,
        'seconds': perf_counter() - start,
        **metadata,
    }
    if hasattr(state, 'algorithm1'):
        result.update(
            schema_version=2,
            algorithm1_status=state.algorithm1.status,
            algorithm1_steps=state.algorithm1.steps,
            obligations=state.obligations,
            equations_by_channel=state.channel_equations,
            coordinate_kind='exact-linear-observable-expressions',
            observable_expressions={name: [
                {'coefficient': str(c), 'gml': rule.getGMLString()}
                for c, rule in expression.values()]
                for name, expression in adapter.expressions.items()},
        )
    if not result['closed']:
        result['status'] = 'UNKNOWN'

    from mfobs.outputs import matrix

    result['matrix'] = (
        [
            [str(value) for value in row]
            for row in matrix(result).tolist()
        ]
        if result['closed']
        else None
    )
    return result
