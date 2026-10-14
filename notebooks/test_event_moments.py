"""Computed pattern events versus a separate finite word-chain reference."""
import mfobs
from mfobs.applications import ABA
from collections import deque

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from mfobs.models.molecular_rop import inputs as paths
import paper_view as view


def word_reference(word, pools, motif, gated, rates, times):
    # No event record or discovered suffix coordinates are consulted here.
    states = [(word, pools['A'], pools['B'], 1)]
    indices = {states[0]: 0}
    pending = deque([0])
    rows, cols, values = [], [], []
    while pending:
        j = pending.popleft()
        w, a, b, opened = states[j]
        transitions = []
        if opened and a:
            transitions.append(((w+'A', a-1, b, opened), rates['k_A']*a))
        if opened and b:
            transitions.append(((w+'B', a, b-1, opened), rates['k_B']*b))
        if gated:
            transitions.append(((w, a, b, 1-opened), rates['a'] if opened else rates['b']))
        rows.append(j); cols.append(j); values.append(-sum(rate for _, rate in transitions))
        for target, rate in transitions:
            if target not in indices:
                indices[target] = len(states); pending.append(len(states)); states.append(target)
            rows.append(indices[target]); cols.append(j); values.append(rate)
    matrix = coo_matrix((values, (rows, cols)), shape=(len(states), len(states))).tocsr()
    p0 = np.zeros(len(states)); p0[0] = 1.
    p = expm_multiply(matrix, p0, start=0, stop=times[-1], num=len(times))
    q = np.array([sum(w[i:i+len(motif)] == motif for i in range(len(w)-len(motif)+1))
                  for w, _, _, _ in states])
    mean = p @ q
    return mean, p @ (q*q)-mean**2


@pytest.mark.parametrize('motif,gated,word', [
    ('ABA', True, 'ABAB'), ('ABBA', True, 'ABAB'),
    ('AAAA', False, 'AAA'), ('BBA', False, 'ABAB'),
])
def test_changed_motif_matches_independent_word_chain(motif, gated, word):
    backend1 = ABA(tuple((''.join(word) for word in paths.configured_variant(motifs=(motif,), gated=gated).motifs)) or paths.configured_variant(motifs=(motif,), gated=gated).alphabet, alphabet=paths.configured_variant(motifs=(motif,), gated=gated).alphabet, gated=paths.configured_variant(motifs=(motif,), gated=gated).gated)
    record = mfobs.analyze(backend1.specification(), backend=backend1, linear_limit=8).record()
    pools = {'A': 3, 'B': 2}
    initial = view.word_coordinates(record, word, pools)
    rates = dict(beta_A=0., beta_B=0., k_A=.8, k_B=.7)
    if gated:
        rates.update(a=.4, b=.6)
    times = np.linspace(0, 4, 21)
    solution = view.integrate_events(record, initial, rates, times, query='C_'+motif)
    mean, variance = word_reference(word, pools, motif, gated, rates, times)
    np.testing.assert_allclose(solution.mean, mean, atol=2e-11)
    np.testing.assert_allclose(solution.variance, variance, atol=2e-11)
    assert solution.diagnostics['exit_probability_bound'] == 0.


def test_rate_sliders_preserve_the_initial_plot_settings(monkeypatch):
    import visual_view as visual
    seen = {}
    def curve(record, initial, rates, **options):
        seen.update(rates=rates, **options)
        return 'plot'
    monkeypatch.setattr(visual, 'aba_curve', curve)
    monkeypatch.setattr(visual, '_live', lambda controls, callback, baseline: callback)
    options = dict(end=100., paths=2000, seed=12345, truncation_tolerance=2e-6, max_states=170000)
    callback = visual.aba_sliders({}, {}, {'k_A': .8, 'k_B': .7, 'beta_A': 2.}, {'options': options})
    assert callback(k_A=.2, k_B=1.8) == 'plot'
    assert {key: seen[key] for key in options} == options
    assert seen['rates'] == {'k_A': .2, 'k_B': 1.8, 'beta_A': 2.}


def test_live_controls_update_visible_plot_report_errors_and_disconnect(monkeypatch):
    import ipywidgets as widgets
    import visual_view as visual
    rendered=[]
    monkeypatch.setattr(visual,'display',lambda *args,**kwargs:None)
    baseline=visual.CurveData(options={'end':20.,'paths':100,'seed':123},
                              rates={'k_A':.9,'k_B':.01,'a':.6,'b':.4})
    baseline.display_handle=type('Handle',(),{'update':lambda self,fig:rendered.append(fig)})()
    slider=widgets.FloatSlider(value=.01,min=0,max=2,step=.01)
    calls=[]
    def compute(**rates):
        calls.append(rates)
        assert slider.disabled
        assert 'Calculating' in panel.children[-1].value
        if rates['k_B']==.03:raise ValueError('budget exhausted <retry>')
        return visual.plt.figure(),{'rates':{**baseline['rates'],**rates}}
    panel=visual._live({'k_B':slider},compute,baseline)
    assert all(s in panel.children[-1].value for s in ['T=20','RUNS=100','a=0.6','b=0.4'])
    slider.value=.02
    assert len(rendered)==1 and baseline['rates']['k_B']==.02
    assert 'Updated' in panel.children[-1].value and not slider.disabled
    slider.value=.03
    assert 'Calculation failed: ValueError: budget exhausted &lt;retry&gt;' in panel.children[-1].value
    assert len(rendered)==1 and baseline['rates']['k_B']==.02 and not slider.disabled
    visual.close_controls(panel)
    slider.value=.04
    assert len(calls)==2 and slider.disabled


def test_aba_plot_recreates_controls_without_overwriting_user_settings():
    from pathlib import Path
    import nbformat
    path=Path(__file__).parent/'examples/02_aba.ipynb'
    n=nbformat.read(path,4)
    plotting=next(c.source for c in n.cells if c.cell_type=='code' and 'controls = visual.aba_sliders' in c.source)
    assert plotting.index('close_controls') < plotting.index('show_curve') < plotting.index('aba_sliders')
    assert 'T = ' not in plotting and 'RUNS = ' not in plotting
    assert 'end=T, paths=RUNS' in plotting
    assert sum('controls = visual.aba_sliders' in c.source for c in n.cells)==1


def test_public_aba_notebook_equations_and_cutoff_check(monkeypatch):
    """Execute the actual consumer cells: motif labels are not event coordinates."""
    import json
    from pathlib import Path
    from mfobs import analyze, ABA
    import visual_view as visual

    notebook = json.loads((Path(__file__).parent / 'examples/02_aba.ipynb').read_text(encoding='utf-8'))
    sources = [''.join(cell['source']) for cell in notebook['cells'] if cell['cell_type'] == 'code']
    backend = ABA.motif('ABA', gated=True)
    record = analyze(backend.specification(), backend=backend, linear_limit=8).record()
    initial = view.word_coordinates(record, 'ABAB', {'A': 2, 'B': 1}, opened=True)
    rates = dict(beta_A=0., beta_B=0., k_A=.9, k_B=.7, a=.4, b=.6)
    figure, curves = visual.aba_curve(record, initial, rates, end=.2, paths=8)
    visual.plt.close(figure)
    monkeypatch.setattr(view, 'display', lambda *args: None)
    monkeypatch.setattr(visual, 'display', lambda *args: None)
    namespace = dict(
        result=record, linear=record['linear'], MOTIF='ABA', initial=initial, curves=curves,
        view=view, visual=visual, display=lambda *args: None,
        NOTEBOOKS=Path(__file__).parent,
    )
    equations = next(source for source in sources if 'view.show_equations(' in source)
    exec(compile(equations, 'ABA equations cell', 'exec'), namespace)
    # The expensive numerical refinement belongs in this regression, not in every reading.
    reference = view.integrate_events(record, initial, rates, curves['time'],
                                      query=record['support_discovery']['selected_queries'][0], truncation_tolerance=1e-8)
    for moment, error in [('mean', 'mean_error_bound'), ('variance', 'variance_error_bound')]:
        difference = np.max(np.abs(curves[moment] - getattr(reference, moment)))
        allowance = curves['numerical'][error] + reference.diagnostics[error] + 1e-9
        assert difference <= allowance


def test_nac_and_remy_controls_preserve_nondefault_settings(monkeypatch):
    import visual_view as visual
    def live(controls, callback, baseline):
        return callback(**{key: widget.value for key, widget in controls.items()})
    monkeypatch.setattr(visual, '_live', live)
    monkeypatch.setattr(visual, 'context_nac_curve', lambda results, **kw: kw)
    monkeypatch.setattr(visual, 'remy_curve', lambda result, **kw: kw)
    assert visual.context_nac_sliders({}, {'initial': {}}, mask=6, rate=3., end=8.) == {'mask':6, 'initial': {}, 'rate':3., 'end':8.}
    assert visual.remy_sliders({}, {}, rate=.8, end=9.) == {'rate':.8, 'end':9.}
    assert visual.remy_sliders({}, {}, rate=0., end=.1) == {'rate':0., 'end':.1}
