"""Fresh ABA records retain their input rules and exact joint jumps."""
from copy import deepcopy
from hashlib import sha256

import pytest
from mfobs import analyze, ABA
from mfobs.models.rop_exploration import expand
CASES = (("aba", False), ("aba-gated", True))



@pytest.fixture(scope='module')
def records(aba_records):
    return {name: {"analysis": result} for name, result in aba_records.items()}


@pytest.mark.parametrize('identifier,gated', CASES)
def test_inputs_and_events_belong_to_one_model(records, identifier, gated):
    record = records[identifier]
    result = record['analysis']
    support = result['support_discovery']
    assert support['selected_queries'] == ['C_ABA']
    assert len(support['patterns']) == 5 + gated
    assert result['pjc']['event_class_count'] == 8 + 2*gated
    assert result['linear']['status'] == 'UNKNOWN_DIMENSION_LIMIT'
    assert result['pjc']['status'] == 'CLOSED'
    actual_inputs = ABA.motif('ABA', gated=gated).rule_gmls()
    for rule in result['contract']['normalized_rules']:
        gml = actual_inputs[rule['id']]
        assert support['source_rules'][rule['id']]['gml'] == gml
        assert rule['gml_sha256'] == sha256(gml.encode()).hexdigest()


def test_literal_application_counts_both_free_tokens_and_checks_the_joint_jump(records):
    expansion = expand(records['aba']['analysis'], depth=1)
    before = expansion['states'][0]
    event = next(t for t in expansion['transitions'] if t['rule'] == 'insert_A')
    after = expansion['states'][event['target']]
    assert (before['word'], after['word'], event['multiplicity']) == ('ABAB', 'ABABA', 2)
    assert [after['coordinates'][n]-v for n,v in before['coordinates'].items()] == [1,-1,0,-1,1]
    damaged = deepcopy(records['aba']['analysis'])
    term = next(t for t in damaged['pjc']['operator']['terms'] if t['delta'].get('C_ABA') == 1)
    term['delta']['C_ABA'] = 2
    with pytest.raises(AssertionError):
        expand(damaged, depth=1)


def test_closed_gate_blocks_append_but_allows_feed_and_opening(records):
    result = expand(records['aba-gated']['analysis'], opened=False, depth=1)
    assert {t['rule'] for t in result['transitions']} == {'feed_A', 'feed_B', 'donor_unbind'}


def test_changed_query_uses_literal_graph_counts_without_a_preset():
    backend = ABA.motif('BABBA')
    result = analyze(backend.specification(), backend=backend, linear_limit=24).record()
    expansion = expand(result, word='BABB', depth=1)
    event = next(t for t in expansion['transitions'] if t['rule'] == 'insert_A')
    assert expansion['states'][event['target']]['coordinates']['C_BABBA'] == 1
    assert result['pjc']['event_class_count'] == 12


