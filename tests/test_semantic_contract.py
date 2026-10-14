from types import SimpleNamespace

import pytest

from mfobs.core.semantic_contract import (
    RulePreservation,
    SemanticContract,
    theorem_records_from_input,
)
from mfobs.core.synthesis import Reduction
from mfobs.core.worklist import _semantic_contract
from mfobs.inputs import InputError


def preservation(rule='reaction', status='VERIFIED'):
    return RulePreservation(
        rule=rule,
        status=status,
        evidence='theorem-argument',
        scope='every graph in the declared test language',
        certificate={'proof_document': 'docs/example.md'},
    )


def test_reduction_status_is_closed_and_unknown_needs_a_reason():
    with pytest.raises(InputError, match='Invalid reduction status'):
        Reduction('approximate', {})
    with pytest.raises(InputError, match='explicit reason'):
        Reduction('unknown', {})


def test_contract_requires_exactly_the_input_rule_set():
    contract = SemanticContract(
        domain='example-language',
        reducer='exact example normal form',
        exact_reduction_scope='all example patterns',
        preservation=(preservation(),),
    )
    with pytest.raises(InputError, match='Preservation rules differ'):
        contract.validate_rules(['reaction', 'second_reaction'])


def test_theorem_input_requires_one_nonempty_claim_per_rule():
    evidence = {
        'level': 'theorem-argument',
        'scope': 'all parameter values',
        'proof_document': 'docs/example.md',
        'rule_claims': {'reaction': 'Preserves the declared family.'},
    }
    records = theorem_records_from_input(['reaction'], evidence)
    assert records[0].certificate['rule_claim'].startswith('Preserves')

    evidence['rule_claims'] = {}
    with pytest.raises(InputError, match='rule claims differ'):
        theorem_records_from_input(['reaction'], evidence)


@pytest.mark.parametrize('case', ['missing', 'unknown'])
def test_required_contract_fails_closed(case):
    adapter = SimpleNamespace(semantic_contract_required=True)
    if case == 'unknown':
        adapter.semantic_contract = SemanticContract(
            domain='example-language',
            reducer='exact example normal form',
            exact_reduction_scope='all example patterns',
            preservation=(preservation(status='UNKNOWN'),),
        )
    state = SimpleNamespace(unknown=[])

    _semantic_contract(adapter, state)

    assert state.unknown
