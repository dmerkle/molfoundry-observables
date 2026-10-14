"""Fresh discovery, changed inputs, and an independent directed-multigraph oracle."""
from mfobs.core.rules import rule_graph
from mfobs.core.graph import make_graph, edge_rows, node_labels
import copy
from fractions import Fraction
import json
from pathlib import Path

import pytest
from mfobs.models.dna_walker.analysis import analyze
from mfobs.outputs import dna_walker_readout
from mfobs.models.dna_walker.inputs import identity
from mfobs.models.dna_walker.counting import reduce_count
from mfobs.checks.dna_walker import validate

from mfobs.core.linear_closure import UnsupportedImage

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT/'notebooks/visual/evidence/dna-walker'

@pytest.fixture(scope='module')
def inputs():
    return json.loads((ROOT/'examples/dna-walker/input.json').read_text(encoding='utf-8'))

@pytest.fixture(scope='module')
def recovered():
    return json.loads((EVIDENCE/'recovered/walker-export.json').read_text(encoding='utf-8'))

@pytest.fixture(scope='module')
def result(inputs):
    result = analyze(inputs['rules'], inputs['query'])
    result.update(dna_walker_readout(result,inputs['rules']))
    return result


def test_input_contains_only_original_rules_and_single_query(inputs,recovered):
    assert set(inputs) == {'rules','query'}
    assert inputs == {'rules':recovered['rules'],'query':recovered['extended_query']}


def test_fresh_discovery_matches_retained_comparison(result):
    archived = json.loads((EVIDENCE/'result.json').read_text(encoding='utf-8'))
    assert result['status'] == 'CLOSED'
    assert list(result['basis']) == ['E','O1']
    assert result['origins']['O1']['origin'] == 'computed-channel-residual'
    assert result['raw_composites'] == 800
    assert result['ordered_composition_calls'] == 16
    assert result['synthesis_state_evaluations'] == 0
    rational = lambda rows:{k:{n:Fraction(c) for n,c in row.items()} for k,row in rows.items()}
    assert rational(result['raw_images']) == rational(archived['raw_images'])
    assert result['comparison']['physical_equations'] == archived['comparison']['physical_equations']
    assert result['comparison']['stationary_velocity'] == archived['comparison']['stationary_velocity']


def test_independent_literal_applications(result,inputs,recovered):
    checked = validate(result,inputs['rules'],{'E':recovered['extended_query'],'Craw':recovered['compressed_query']})
    assert checked['hosts'] == 136
    assert checked['literal_multigraph_rule_applications'] == 408
    assert checked['generator_images_checked'] == 1088
    assert checked['original_vs_incidence_observable_counts'] == 272
    assert checked['individual_projected_pattern_checks'] > 0


@pytest.mark.parametrize('selected',['extended_query','compressed_query'])
@pytest.mark.parametrize('scale',[1,-3])
def test_changed_query_and_scale_keep_physical_result(selected,scale,recovered,result):
    changed = analyze(recovered['rules'],recovered[selected],query_name='selected',query_scale=scale)
    assert changed['status'] == 'CLOSED'
    assert list(changed['basis']) == ['selected','O1']
    assert list(changed['basis']['selected'].values()) == [str(scale)]
    readout = dna_walker_readout(changed,recovered['rules'])
    assert readout['raw_images'] == result['raw_images']


def test_remove_extension_rules_changes_discovery(inputs):
    rules = [r for r in inputs['rules'] if r['name'] in ('kFC','kBC')]
    result = analyze(rules,inputs['query'])
    assert result['status'] == 'CLOSED'
    assert list(result['basis']) == ['E']  # No preloaded compressed coordinate.


def test_renamed_channels_are_not_interpreted(inputs,result):
    rules = copy.deepcopy(inputs['rules'])
    for i,r in enumerate(rules): r['name'] = 'channel_'+str(i)
    fresh = analyze(rules,inputs['query'])
    assert fresh['basis'] == result['basis']
    assert [s['coordinates'] for s in fresh['steps']] == [s['coordinates'] for s in result['steps']]


def test_budget_exit_does_not_claim_closure(inputs):
    result = analyze(inputs['rules'],inputs['query'],maximum_dimension=1)
    assert result['status'] == 'UNKNOWN_DIMENSION_LIMIT'
    assert result['obligations'] and result['pending']


def test_changing_scaffold_is_outside_declared_language(inputs):
    rule = copy.deepcopy(inputs['rules'][0])
    rule['preserved_edges'] = [pair for pair in rule['preserved_edges'] if pair[0] != 4]
    with pytest.raises(UnsupportedImage,match='preserve the scaffold'):
        analyze([rule],inputs['query'])


def test_incomplete_context_is_unsupported_not_zero(inputs):
    graph = rule_graph(identity(inputs['query'],'q'))
    edge = edge_rows(graph)[0]
    with pytest.raises(UnsupportedImage,match='Incomplete edge'):
        reduce_count(make_graph(node_labels(graph),[e for e in edge_rows(graph) if e != edge]))


def test_longer_scaffold_context_is_not_certified_by_finite_samples(inputs):
    query = copy.deepcopy(inputs['query'])
    for i in range(3,6):
        query['nodes'].append({'id':'c'+str(i),'label':'chain'})
        query['edges'].append({'id':10+i,'source':'c'+str(i-1),'target':'c'+str(i)})
    with pytest.raises(UnsupportedImage,match='shorter than five'):
        reduce_count(rule_graph(identity(query,'long')))


def test_additional_rule_conditions_are_not_silently_dropped(inputs):
    rules = copy.deepcopy(inputs['rules'])
    rules[0]['no_edges'] = [['b','c2']]
    with pytest.raises(UnsupportedImage,match='application conditions'):
        analyze(rules,inputs['query'])
