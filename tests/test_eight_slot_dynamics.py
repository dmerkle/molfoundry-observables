"""Regression checks of eight-site counts, joint jumps and finite dynamics."""
from mfobs.core.graph import make_graph
import json
from pathlib import Path
import sys

import numpy as np
import networkx as nx
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'tests/data'
sys.path.insert(0, str(ROOT/'scripts'))
from mfobs.checks.eight_slots import TOPOLOGIES
from mfobs.core.rules import is_identity
from conftest import regression_record
import molfoundry as mf

from mfobs.core.mixtures import fast_count


@pytest.fixture(scope='module')
def report():
    return regression_record('eight-dynamics')


def test_rule_spans_unchanged_and_observables_are_strict_identities(report):
    assert len(report['patterns']) == 48
    for p in report['patterns']:
        for term in p['terms']:
            gml=(OUT/'eight-observables'/Path(term['gml']).name).read_text(encoding='utf-8')
            assert is_identity(mf.Rule.fromGMLString(gml, add=False))


@pytest.mark.parametrize('tag', list(TOPOLOGIES))
def test_all_initial_coordinates_recounted_from_actual_graphs(report, tag):
    from mfobs.models.eight_slots_common import _representative
    host = _representative(tuple((0,) * n for n in TOPOLOGIES[tag]))
    initial = report['initial_conditions'][tag]['observable_initial_values']
    from fractions import Fraction
    for p in report['patterns']:
        value=Fraction()
        for term in p['terms']:
            rule=mf.Rule.fromGMLString((OUT/'eight-observables'/Path(term['gml']).name).read_text(encoding="utf-8"), add=False)
            value+=Fraction(term['coefficient'])*fast_count(rule,host)
        assert value==Fraction(str(initial[p['id']]))
    assert [Fraction(str(initial[n])) for n in ('P0','P1','P2')]==[0,8,0]


def test_depth_three_public_expansion_and_disconnected_components(report):
    direct = report['direct_expansion']
    assert direct['public_apply_calls'] == direct['independent_target_set_checks'] > 0
    graphs = {g['id']: g for g in direct['graphs']}
    assert [len(c['states']) for c in direct['cases']] == [11, 8, 13]
    for case in direct['cases']:
        levels = {s['graph']: s['depth'] for s in case['states']}
        assert max(levels.values()) == 3 and levels[case['initial_graph']] == 0
        for gid in levels:
            g = graphs[gid]; literal = make_graph(g['graph']['nodes'], g['graph']['edges'])
            assert tuple(sorted(map(len, nx.connected_components(nx.Graph([(e.source.id, e.target.id) for e in literal.edges]))))) == tuple(3*n for n in TOPOLOGIES[case['id']])
        for edge in case['edges']:
            assert levels[edge['target']] <= levels[edge['source']]+1
            assert levels[edge['source']] < 3 and edge['multiplicity'] > 0
            change = graphs[edge['target']]['occupied'] - graphs[edge['source']]['occupied']
            assert change == (-1 if edge['rule'] == 'unbind' else 1)
        # Two disconnected identical empty components must not be collapsed to
        # one: all starts have eight basal matches, independent of topology.
        assert sum(e['multiplicity'] for e in case['edges'] if e['source'] == case['initial_graph'] and e['rule'] == 'basal_bind') == 8


def test_saved_closure_and_first_composition_counters(report):
    assert report['synthesis']['closed'] and not report['synthesis']['pending'] and not report['synthesis']['unknown']
    assert [int(row['after']) for row in report['discovery_rounds']][-2:] == [48, 48]
    assert all(int(row['after']) >= int(row['before']) for row in report['discovery_rounds'])
    assert [(row['left_raw'], row['right_raw']) for row in report['first_compositions']] == [(1, 2), (2, 3), (2, 1)]
    assert report['verification']['synthesis_used_no_host_or_orbit_evaluations']


def test_trajectory_replay_count_identity_and_topology_difference(report):
    assert report['replayed_occupancy_max_error'] < 1e-9
    dependency = report['verification']['verified_additional_count_dependencies'][0]
    from fractions import Fraction
    relation = np.array([float(Fraction(dependency.get(name, 0))) for name in report['basis']])
    # An exact nullspace vector has arbitrary scale. Test numerical error for
    # a unit-norm relation, independently of that elimination scale.
    relation /= np.linalg.norm(relation, ord=1)
    arrays = {}
    for curve in report['curves']:
        y = np.array([row['values'] for row in curve['rows']])
        assert y.shape == (151, 48) and np.isfinite(y).all()
        assert np.max(abs(y[:, 0]+y[:, 1]-8)) < 1e-9
        assert np.max(abs(y @ relation)) < 1e-8
        assert np.min(y[:, 0]) > -1e-9 and np.max(y[:, 0]) < 8+1e-9
        arrays[curve['case']] = y
    assert max(abs(arrays['8'][:, 0]-arrays['4+4'][:, 0])) > .8
    assert max(abs(arrays['8'][:, 0]-arrays['3+5'][:, 0])) > .8


def test_trajectory_check_accepts_public_record_without_export_metadata(tmp_path):
    from check_eight_slots import numeric_trajectories
    from mfobs.checks import eight_slots as ref
    families = {}
    for tag, lengths in ref.TOPOLOGIES.items():
        states = list(ref.assignments(lengths))
        families[tag] = {'states': states, 'edges': ref.topology_edges(lengths),
                         'counts': {bits: {'P0': len(lengths)} for bits in states}}
    # Component count is invariant under every concrete occupancy transition.
    result = {'basis': ['P0'], 'equations_by_reaction': {'P0': {}}}
    report = numeric_trajectories(result, families, tmp_path, write_observables=False)
    assert report['max_ODE_vs_independent_CTMC_error_all_observables'] < 1e-8
    assert (tmp_path / 'topology-trajectories.svg').is_file()
    assert not (tmp_path / 'observables.md').exists()
