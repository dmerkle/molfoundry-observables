"""Numerical and presentation contracts for the NAC comparison."""
from pathlib import Path
import numpy as np
import pytest
from mfobs.numerics.context_nac import initial_counts, forward_system, trajectory

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'experiments/chlorination/context-nac'


@pytest.fixture(params=['guarded','ablated'])
def record(request, chlorination_records):
    return request.param, chlorination_records[request.param]


def test_generated_event_forward_equations_and_limits(record):
    variant, result = record
    for mask in (0,1,9):
        initial = initial_counts(result, mask)
        solution = trajectory(result, initial, times=[0.,1.,100.])
        assert np.allclose(solution['probabilities'].sum(1), 1)
        assert solution['Z'][0] == mask.bit_count()
        assert solution['Z'][-1] == pytest.approx(2 if variant=='guarded' else 4, abs=1e-8)
        assert solution['variance'][-1] == pytest.approx(0, abs=1e-8)
        assert solution['blocked'][-1] == pytest.approx(1 if variant=='guarded' else 0, abs=1e-8)
        stopped = trajectory(result, initial, times=[0.,3.], rate=0)
        assert np.all(stopped['Z'] == initial['Z'])
        assert np.all(stopped['variance'] == 0)


def test_forward_search_refuses_truncated_results(record):
    _, result = record
    with pytest.raises(RuntimeError, match='budget exhausted'):
        forward_system(result, initial_counts(result), maximum_states=1)


def test_slider_reuses_initial_counts_without_calling_graph_engine(monkeypatch, chlorination_records):
    import importlib.util
    import matplotlib.pyplot as plt
    import mfobs.numerics.context_nac as numerical
    monkeypatch.syspath_prepend(str(ROOT/'notebooks'))
    spec = importlib.util.spec_from_file_location('nac_visual_view', ROOT/'notebooks/visual_view.py')
    visual = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(visual)
    records = chlorination_records
    fig, baseline = visual.context_nac_curve(records)
    plt.close(fig)
    def no_graph_calls(*args, **kwargs):
        raise AssertionError('Slider updates must only evaluate the numerical model')
    monkeypatch.setattr(numerical, 'initial_counts', no_graph_calls)
    monkeypatch.setattr(visual, '_live', lambda controls, callback, baseline: callback)
    callback = visual.context_nac_sliders(records, baseline)
    fig, changed = callback(rate=0., end=3.)
    plt.close(fig)
    assert np.all(changed['guarded_Z'] == 0)
    assert np.all(changed['ablated_Z'] == 0)

