"""Plot consumers use computed coefficients independently of event extraction."""
from pathlib import Path
import numpy as np
import pytest
from mfobs.numerics.linear import linear_matrix, trajectory
from mfobs.numerics.events import affine_generator

ROOT = Path(__file__).resolve().parents[1]


def test_remy_mean_ode_needs_no_event_partition_and_matches_known_edge_growth(remy_record):
    from copy import deepcopy
    result = deepcopy(remy_record)
    names, a = linear_matrix(result['linear'], {'G': .25})
    enames, e, offset = affine_generator(result, {'rho': .25})
    order = [enames.index(n) for n in names]
    assert np.allclose(a, e[np.ix_(order, order)]) and np.allclose(offset, 0)
    initial = {n: float(n == 'E') for n in names}
    result.pop('pjc')
    curve = trajectory(result['linear'], initial, {'G': .25}, end=3.)
    assert np.allclose(curve['values'][:, names.index('E')], np.exp(curve['time']), rtol=1e-9)


def test_linear_consumer_refuses_incomplete_records_and_missing_channels():
    with pytest.raises(ValueError, match='closed'):
        linear_matrix({'status':'UNKNOWN'}, {})
    r = {'status':'CLOSED','basis':['X'], 'steps':[{'source':'X','channel':'G','coordinates':{'X':'2'}}]}
    with pytest.raises(ValueError, match='channel rates'):
        linear_matrix(r, {})
    with pytest.raises(ValueError, match='nonnegative'):
        linear_matrix(r, {'G': -1})
    r['basis'].append('Y')
    with pytest.raises(ValueError, match='Incomplete'):
        linear_matrix(r, {'G': 1})
