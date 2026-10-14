"""Regression tests for configurations that previously returned false CLOSED."""
from dataclasses import replace
from pathlib import Path
import copy
import json
import shutil

import pytest

from mfobs.core.forbidden import ForbiddenPatternSpecification, IncidentBound
from mfobs.core.worklist import run
from mfobs.inputs import InputError, load, with_observable_queries
from mfobs.query_input import from_arguments
from mfobs.invariants.chain import ChainLanguage
from mfobs.invariants.contracts import validate_language, _profiles, _check_profile

ROOT = Path(__file__).resolve().parents[1]


def gated_path():
    return load(ROOT / 'examples/gated-path/model.json')


@pytest.mark.parametrize('adapter', ['mfobs.models.chain_mixture:Adapter',
                                   'mfobs.core.identity:Adapter'])
def test_connected_overlap_false_closed_is_rejected_before_composition(adapter):
    inp = gated_path()
    inp.data['adapter'] = adapter
    inp.data['composition']['connected'] = True
    with pytest.raises(InputError, match='UNSOUND_OVERLAP_POLICY'):
        run(inp)


def test_connected_policy_also_rejected_at_manifest_load(tmp_path):
    dst = tmp_path / 'example'
    shutil.copytree(ROOT / 'examples/gated-path', dst)
    path = dst / 'model.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    data['composition']['connected'] = True
    path.write_text(json.dumps(data), encoding='utf-8')
    with pytest.raises(InputError, match='UNSOUND_OVERLAP_POLICY'):
        load(path)


def test_incompatible_valence_false_closed_is_rejected_before_composition():
    inp = gated_path()
    original = inp.forbidden_spec()
    incompatible = replace(original, max_valence={**original.max_valence, 'A': 1})
    inp.forbidden_spec = lambda: incompatible
    with pytest.raises(InputError, match='LANGUAGE_FORBIDDEN_MISMATCH.*valence 2 > 1'):
        run(inp)


def test_ring_connected_exception_is_still_digest_bound():
    inp = load(ROOT / 'examples/ring8/model.json')
    assert inp.data['composition']['connected']
    inp.adapter()
    inp.data['rules'][0]['per_match'] = '2'
    with pytest.raises(InputError, match='MODEL_CONTRACT_MISMATCH'):
        inp.adapter()


def test_ring_backend_binding_allows_a_separately_hashed_query_selection():
    inp = load(ROOT / 'examples/ring8/model.json')
    selection = from_arguments((
        f"P0={ROOT / 'examples/ring8/observables/occupied.gml'}",
    ))
    overridden = with_observable_queries(inp, selection)
    assert overridden.backend_digest() == inp.backend_digest()
    assert overridden.definition_digest() != inp.definition_digest()
    overridden.adapter()


@pytest.mark.parametrize('value', [-1, 1.5, '1', True, None, 2**53])
@pytest.mark.parametrize('field', ['bond_orders', 'max_valence', 'incident_bound'])
def test_malformed_numerical_pruning_bounds_are_rejected(field, value):
    spec = gated_path().forbidden_spec()
    with pytest.raises(ValueError, match='nonnegative safe integer'):
        if field == 'bond_orders':
            replace(spec, bond_orders={**spec.bond_orders, '-': value})
        elif field == 'max_valence':
            replace(spec, max_valence={**spec.max_valence, 'A': value})
        else:
            IncidentBound('A', '-', 'A', value, 'neighbour-count')


@pytest.mark.parametrize('change', [
    {'atom_labels': ('A', 'A')},
    {'strict_signature': 'true'},
    {'bond_orders': {'=': 2}},
    {'allowed_bond_pairs': {'-': (('A', 'F'),)}},
    {'allowed_bond_pairs': {'-': (('A',),)}},
    {'max_valence': {'missing-label': 1}},
    {'incident_bounds': (IncidentBound('A', 'missing-edge', 'A', 1, 'bad'),)},
    {'incident_bounds': (IncidentBound('F', '-', 'A', 1, 'bad'),)},
])
def test_signature_inconsistencies_are_rejected(change):
    with pytest.raises(ValueError):
        replace(gated_path().forbidden_spec(), **change)


def test_pair_orientation_certificate_matches_executable_payload():
    spec = ForbiddenPatternSpecification('test', ('A', 'B'), {'-': 1},
        {'-': (('B', 'A'),)}, {'A': 1, 'B': 1})
    assert spec.allowed_bond_pairs['-'] == (('A', 'B'),)
    assert {'kind': 'invalid-bond-signature', 'nodes': ['A', 'B'], 'edge': '-'} not in spec.invalid_bond_patterns()


def test_zero_order_relation_is_supported_with_explicit_incidence_bound():
    spec = ForbiddenPatternSpecification('relation', ('A',), {'-': 0},
        {'-': (('A', 'A'),)}, {'A': 0},
        (IncidentBound('A', '-', 'A', 2, 'neighbours'),))
    assert not spec.minimal_overvalence_patterns()
    assert spec.certificate()['counts']['incident_bound_patterns'] == 1


def path_language():
    tile = {'nodes': [[0, 'C']], 'edges': []}
    return {'chain': {'kind': 'one-port-decorated-path', 'join_bond': '-',
        'zero': tile,
        'left': {**tile, 'incoming': None, 'outgoing': 0},
        'repeat': {**tile, 'incoming': 0, 'outgoing': 0},
        'right': {**tile, 'incoming': 0, 'outgoing': None}},
        'open': tile, 'closed': tile, 'monomer': tile}


def path_theory(maximum):
    return ForbiddenPatternSpecification('path-stars', ('C',), {'-': 1},
        {'-': (('C', 'C'),)}, {'C': maximum})


def test_incoming_and_outgoing_same_port_both_contribute():
    language = path_language()
    # Each isolated tile satisfies valence<=1. A repeat tile needs both joins.
    for tile in ['zero', 'left', 'repeat', 'right']:
        assert language['chain'][tile]['edges'] == []
    with pytest.raises(InputError, match=r'chain.left\|repeat\|right.*valence 2 > 1'):
        validate_language(language, path_theory(1))


def test_local_contexts_include_arbitrarily_far_interior_without_length_sampling(monkeypatch):
    language = path_language()
    theory = path_theory(2)
    original = ChainLanguage.instantiate
    def no_sampling(*args, **kwargs):
        raise AssertionError('Contract must derive local contexts, not instantiate a length cutoff')
    monkeypatch.setattr(ChainLanguage, 'instantiate', no_sampling)
    certificate = validate_language(language, theory)
    assert certificate['status'] == 'VERIFIED_ALL_SIZES'
    assert len(certificate['chain_contexts']) == 8
    assert {'previous': 'repeat', 'tile': 'repeat', 'next': 'repeat'} in certificate['chain_contexts']
    middle = [p for p in certificate['profiles'] if p['context'] == 'chain.repeat|repeat|repeat']
    assert len(middle) == 1 and middle[0]['valence'] == 2
    # A large independent instance corroborates the classification; it is NOT
    # part of the membership proof, and the contract above cannot inspect it.
    monkeypatch.setattr(ChainLanguage, 'instantiate', original)
    host, _ = ChainLanguage(language['chain']).instantiate(1001)
    for profile in _profiles(host, 'independent-long-chain'):
        _check_profile(profile, theory)


@pytest.mark.parametrize('part', ['open', 'closed', 'monomer'])
def test_every_finite_component_is_checked(part):
    language = path_language()
    language[part] = {'nodes': [[0, 'O']], 'edges': []}
    with pytest.raises(InputError, match=f'LANGUAGE_FORBIDDEN_MISMATCH: {part}'):
        validate_language(language, path_theory(2))


def test_incidence_and_join_bond_constraints_are_checked():
    theory = path_theory(2)
    theory = replace(theory, incident_bounds=(IncidentBound('C', '-', 'C', 1, 'at-most-one'),))
    with pytest.raises(InputError, match='incident bound at-most-one: 2 > 1'):
        validate_language(path_language(), theory)
    theory = replace(path_theory(2), allowed_bond_pairs={'-': ()})
    with pytest.raises(InputError, match='forbidden bond'):
        validate_language(path_language(), theory)


@pytest.mark.parametrize('mutation', ['duplicate-node', 'duplicate-edge', 'missing-port', 'bool-port'])
def test_malformed_grammar_fails_before_normalization(mutation):
    language = copy.deepcopy(path_language())
    if mutation == 'duplicate-node':
        language['chain']['repeat']['nodes'] *= 2
    elif mutation == 'duplicate-edge':
        language['monomer'] = {'nodes': [[0, 'C'], [1, 'C']], 'edges': [[0, 1, '-'], [1, 0, '-']]}
    elif mutation == 'missing-port':
        del language['chain']['repeat']['incoming']
    else:
        language['chain']['repeat']['incoming'] = False
    with pytest.raises(InputError, match='LANGUAGE_CONTRACT'):
        validate_language(language, path_theory(2))
