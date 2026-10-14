"""Uniform evidence contract between invariant reducers and the core.

The core does not know chemistry, cycles, chains, or finite state spaces.  It
does need an inspectable answer to two questions before reporting a conditional
closure: what reductions are exact, and why do the rules preserve the declared
state language?  These small immutable records are that boundary.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

from mfobs.inputs import InputError


PRESERVATION_STATUSES = {'VERIFIED', 'UNKNOWN', 'VIOLATION'}
EVIDENCE_LEVELS = {'machine-check', 'exhaustive-finite', 'theorem-argument'}


@dataclass(frozen=True)
class RulePreservation:
    rule: str
    status: str
    evidence: str
    scope: str
    certificate: dict

    def __post_init__(self):
        if self.status not in PRESERVATION_STATUSES:
            raise InputError(f'Invalid preservation status for {self.rule}')
        if self.evidence not in EVIDENCE_LEVELS:
            raise InputError(f'Invalid preservation evidence for {self.rule}')
        if (
            not self.rule
            or not self.scope
            or not isinstance(self.certificate, dict)
        ):
            raise InputError(f'Malformed preservation record for {self.rule}')

    def record(self):
        return asdict(self)


@dataclass(frozen=True)
class SemanticContract:
    """Reducer scope and one preservation judgement per input rule."""

    domain: str
    reducer: str
    exact_reduction_scope: str
    preservation: tuple[RulePreservation, ...]

    def __post_init__(self):
        if not self.domain or not self.reducer or not self.exact_reduction_scope:
            raise InputError('Semantic contract needs domain, reducer and scope')
        rules = [item.rule for item in self.preservation]
        if not rules or len(rules) != len(set(rules)):
            raise InputError(
                'Semantic contract needs one unique preservation record per rule'
            )

    def validate_rules(self, expected):
        expected = list(expected)
        actual = [item.rule for item in self.preservation]
        if set(actual) != set(expected) or len(actual) != len(expected):
            raise InputError(
                f'Preservation rules differ: expected={expected}, actual={actual}'
            )

    @property
    def allows_conditional_closure(self):
        return all(item.status == 'VERIFIED' for item in self.preservation)

    def record(self):
        levels = sorted({item.evidence for item in self.preservation})
        return {
            'schema': 'mfobs-semantic-contract-v1',
            'domain': self.domain,
            'reducer': self.reducer,
            'reduction_outcomes': ['EXACT', 'ZERO', 'UNKNOWN'],
            'exact_reduction_scope': self.exact_reduction_scope,
            'preservation': {
                item.rule: item.record() for item in self.preservation
            },
            'preservation_evidence_levels': levels,
            'all_rules_verified': self.allows_conditional_closure,
            'allows_conditional_closure': self.allows_conditional_closure,
            'claim': (
                'Closure is conditional on the declared state language and '
                'the listed preservation evidence; it is not an unconditional '
                'claim over all finite graphs.'
            ),
        }


def machine_preservation(rule, scope, certificate, passed=True):
    status = 'VERIFIED' if passed else 'VIOLATION'
    return RulePreservation(rule, status, 'machine-check', scope, certificate)


def theorem_preservation(rule, scope, certificate):
    return RulePreservation(
        rule,
        'VERIFIED',
        'theorem-argument',
        scope,
        certificate,
    )


def exhaustive_preservation(rule, scope, certificate, passed=True):
    status = 'VERIFIED' if passed else 'VIOLATION'
    return RulePreservation(
        rule,
        status,
        'exhaustive-finite',
        scope,
        certificate,
    )


def theorem_records_from_input(
    rule_ids: Iterable[str],
    evidence: dict,
    shared_certificate: dict | None = None,
):
    """Validate explicit human-proof metadata and create per-rule records.

    The function does not pretend to check the mathematical proof.  It makes
    its status, scope, source document, and rule-specific claims impossible to
    omit or accidentally inherit from a different model using the same
    reducer implementation.
    """
    required = {'level', 'scope', 'proof_document', 'rule_claims'}
    if not isinstance(evidence, dict) or set(evidence) != required:
        raise InputError(
            'preservation_evidence must contain exactly level, scope, '
            'proof_document and rule_claims'
        )
    if evidence['level'] != 'theorem-argument':
        raise InputError('Only theorem-argument input evidence is supported here')
    if not isinstance(evidence['scope'], str) or not evidence['scope']:
        raise InputError('preservation_evidence.scope must be nonempty')
    valid_document = (
        isinstance(evidence['proof_document'], str)
        and bool(evidence['proof_document'])
    )
    if not valid_document:
        raise InputError('preservation_evidence.proof_document must be nonempty')

    expected = list(rule_ids)
    claims = evidence['rule_claims']
    if not isinstance(claims, dict) or set(claims) != set(expected):
        raise InputError(
            f'preservation rule claims differ: expected={expected}, '
            f'actual={sorted(claims) if isinstance(claims, dict) else claims}'
        )

    common = dict(shared_certificate or {})
    common['proof_document'] = evidence['proof_document']
    records = []
    for rule_id in expected:
        claim = claims[rule_id]
        if not isinstance(claim, str) or not claim:
            raise InputError(f'Preservation claim for {rule_id} must be nonempty')
        records.append(
            theorem_preservation(
                rule_id,
                evidence['scope'],
                {**common, 'rule_claim': claim},
            )
        )
    return tuple(records)
