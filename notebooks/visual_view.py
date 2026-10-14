"""Shared plots and controls for the example notebooks.

Algorithms remain in mfobs and examples. This module only records runs,
draws their inputs/results, and calls the existing numerical consumers.
"""
from __future__ import annotations
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timezone
from hashlib import sha256
from html import escape
from pathlib import Path
import importlib.metadata
import json
import platform
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import HTML, display
import paper_view as view

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent

class Run:
    def __init__(self, example):
        from mfobs.provenance import runtime_fingerprint
        self.path = HERE / 'runs' / (example + '-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S-%fZ'))
        self.path.mkdir(parents=True)
        sources = sorted((PACKAGE / 'src/mfobs').rglob('*.py')) + sorted((PACKAGE / 'artifact').rglob('*.py'))
        sources += [Path(__file__), HERE / 'paper_view.py', HERE / 'inspection.py']
        sources.append(HERE / 'visual_suite.py')
        sources += sorted((PACKAGE / 'src/mfobs/data').glob('*.json'))
        versions = {}
        for name in ['molfoundry', 'numpy', 'scipy', 'sympy', 'ipywidgets']:
            versions[name] = importlib.metadata.version(name)
        try:
            revision = subprocess.check_output(['git','-C',str(PACKAGE),'rev-parse','HEAD'], text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            revision = 'source snapshot'
        self.metadata = dict(example=example, source_revision=revision,
            timestamp_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
            versions=versions, runtime=runtime_fingerprint(),
            source_sha256={p.relative_to(PACKAGE).as_posix(): sha256(p.read_bytes()).hexdigest() for p in sources})
        (self.path / 'environment.json').write_text(json.dumps(self.metadata, indent=2)+'\n', encoding='utf-8')

    def fresh_path(self, name, suffix):
        path=self.path / (name + suffix)
        if path.exists():
            path=self.path / (name + '-' + datetime.now(timezone.utc).strftime('%H%M%S-%f') + suffix)
        return path

    @contextmanager
    def capture(self, name):
        # Preserve detailed progress without filling the notebook with engine logs.
        with self.fresh_path(name, '.log').open('x', encoding='utf-8') as stream, redirect_stdout(stream):
            yield



def rules(items, *, identity=False):
    """Display each rule or query with MolFoundry's public SVG renderer."""
    from mfobs.visual_links import label_legend
    items = list(items)
    for title, gml in items:
        view.show_molfoundry_rule(gml, title, identity=identity, legend=False)
    legend = label_legend("\n".join(gml for _, gml in items))
    if legend:
        display(HTML("<p>" + escape(legend) + "</p>"))


def table(rows, index=None):
    frame = pd.DataFrame(rows)
    if index is not None: frame = frame.set_index(index)
    display(frame)


def support(record):
    patterns = record['support_discovery']['patterns']
    selected = set(record['support_discovery']['selected_queries'])
    def title(name,p):
        origin = p['origin']
        source = 'query' if name in selected else ('activity: '+origin.get('rule','') if origin['kind']=='rule-activity'
                  else 'from '+origin.get('pattern','')+' / '+origin.get('rule',''))
        meaning = {'pool':'free tokens','suffix':'path ending','motif':'occurrences','total':'total tokens','gate':'open gate'}[p['kind']]
        return f"{name}: {meaning}\n{source}"
    rules([(title(n,p),p['gml']) for n,p in patterns.items() if n not in selected], identity=True)


def aba_curve(record, initial, rates, *, end=4., paths=256, seed=20260924,
              truncation_tolerance=1e-6, max_states=150000):
    if not np.isfinite(end) or end <= 0 or not isinstance(paths,int) or paths < 2:
        raise ValueError('Use a positive end time and at least two SSA paths')
    times=np.linspace(0,end,121)
    names=view.event_system(record).names
    query=record['support_discovery']['selected_queries'][0]; col=names.index(query)
    from mfobs.numerics.events import sample_moments
    data = sample_moments(record, initial, rates, times, paths=paths, seed=seed,
                          truncation_tolerance=truncation_tolerance, max_states=max_states)
    samples = data['samples']
    fig,axes=plt.subplots(1,3,figsize=(12.5,3.8),layout='constrained')
    for path in samples[:8,:,col]: axes[0].step(times,path,where='post',alpha=.6,lw=1.2)
    values=samples[:,:,col]; mean=values.mean(0); error=1.96*values.std(0,ddof=1)/np.sqrt(paths)
    axes[1].plot(times,data['mean'],color='#9c3b35',lw=2,label='Computed (FSP)')
    axes[1].step(times,mean,where='post',color='#2f6595',lw=1.4,label=f'SSA mean ({paths} paths)')
    axes[1].fill_between(times,mean-error,mean+error,color='#2f6595',alpha=.16,label='±1.96 MC standard errors')
    variance=values.var(0,ddof=1)
    axes[2].plot(times,data['variance'],color='#9c3b35',lw=2,label='Computed (FSP)')
    axes[2].step(times,variance,where='post',color='#2f6595',lw=1.4,label=f'SSA variance ({paths} paths)')
    for ax in axes: ax.set(xlabel='Time');ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel(query)
    for ax,title in zip(axes,('Individual paths','Expectation','Variance')):ax.set_title(title)
    for ax in axes[1:]:ax.legend(fontsize=8)
    fig.suptitle(', '.join(f'{key}={value:g}' for key,value in rates.items() if key.startswith('k_')) + f'; initial {query}={initial[query]}',fontsize=11)
    bound=data['numerical']
    fig.supxlabel(f"FSP truncation error over [0, {end:g}]: mean ≤ {bound['mean_error_bound']:.1e}, "
                  f"variance ≤ {bound['variance_error_bound']:.1e} · {bound['states']:,} states",
                  fontsize=9)
    return fig, data


class CurveData(dict):
    """Numerical baseline plus the notebook's updateable PNG display handle."""


def show_curve(result):
    fig,data=result
    data=CurveData(data)
    data.display_handle=display(fig,display_id=True)
    plt.close(fig)
    return data


def close_controls(panel):
    """Disconnect an older panel before replacing the plot it controlled."""
    if panel is not None:
        dispose = getattr(panel, '_dispose_controls', None)
        if dispose is not None:
            dispose()
        panel.close()


def _live(controls, callback, baseline):
    import ipywidgets as w
    options = baseline.get('options', {})
    setting_text = ', '.join(f'{label}={options[key]}' for key,label in
                            [('end','T'),('paths','RUNS'),('seed','seed')]
                            if key in options)
    fixed_rates = ', '.join(f'{key}={value:g}' for key,value in baseline.get('rates',{}).items()
                            if key not in controls)
    context = escape(' · '.join(text for text in (setting_text, fixed_rates) if text))
    status=w.HTML('<small>Ready. Release a slider to update the plot above. '+context+'</small>')
    def update(change):
        status.value='<b>Calculating…</b><br><small>'+context+'</small>'
        for widget in controls.values(): widget.disabled=True
        fig=None
        try:
            fig,data=callback(**{key:widget.value for key,widget in controls.items()})
            baseline.display_handle.update(fig)
            baseline.update(data)
            changed=', '.join(f'{key}={widget.value:g}' for key,widget in controls.items())
            status.value='<small>Updated: '+escape(changed)+'. '+context+'</small>'
        except Exception as error:
            status.value=('<b>Calculation failed: '+escape(type(error).__name__+': '+str(error))+
                          '</b><br><small>The previous plot is unchanged. '+context+'</small>')
        finally:
            if fig is not None: plt.close(fig)
            for widget in controls.values(): widget.disabled=False
    for widget in controls.values(): widget.observe(update,names='value')
    panel=w.VBox([w.HBox(list(controls.values())),status])
    def dispose():
        for widget in controls.values():
            widget.unobserve(update,names='value')
            widget.disabled=True
        status.value='<small>These controls have been replaced by the new settings cell.</small>'
    panel._dispose_controls=dispose
    bundle=panel._repr_mimebundle_()
    bundle['text/plain']='Run All in JupyterLab to activate the sliders.'
    bundle['text/html']='<p><em>Run All in JupyterLab to activate the sliders.</em></p>'
    display(bundle,raw=True)
    return panel


def aba_sliders(record,initial,rates,baseline):
    import ipywidgets as w
    controls={key:w.FloatSlider(description=key,value=rates[key],min=0.,max=max(2.,rates[key]),step=.01,
        continuous_update=False,readout_format='.2f',layout=w.Layout(width='45%')) for key in rates if key.startswith('k_')}
    # Inherit settings from the plot itself; a slider must not restore defaults.
    options=dict(baseline['options'])
    return _live(controls,lambda **values:aba_curve(record,initial,{**rates,**values},**options),baseline)


def remy_curve(result, *, rate=.25, end=3.):
    from mfobs.numerics.linear import trajectory
    initial = {name: float(name == 'E') for name in result['linear']['basis']}
    curves = trajectory(result['linear'], initial, {'G': rate}, end=end)
    fig, ax = plt.subplots(figsize=(10, 3.5), layout='constrained')
    for i, name in enumerate(curves['names']):
        ax.plot(curves['time'], curves['values'][:, i], label=name)
    ax.set(xlabel='Time', ylabel='Expected count', title=f'Closed expectation equations · κ={rate:g}')
    ax.legend(); ax.spines[['top', 'right']].set_visible(False)
    return fig, curves


def remy_sliders(result, baseline, *, rate=.25, end=3.):
    import ipywidgets as widgets
    controls = {'rate': widgets.FloatSlider(value=rate, min=0, max=max(.5, rate), step=.05, description='κ', continuous_update=False),
                'end': widgets.FloatSlider(value=end, min=min(.5, end), max=max(5, end), step=.5, description='Time', continuous_update=False)}
    return _live(controls, lambda **values: remy_curve(result, **values), baseline)


def context_nac_curve(results, *, mask=0, rate=1., end=3., initial=None):
    from mfobs.numerics.context_nac import initial_counts, trajectory
    times = np.linspace(0, end, 201)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5), layout='constrained')
    if initial is None:
        initial = {variant: initial_counts(result, mask) for variant, result in results.items()}
    data = {'time': times, 'initial': initial}
    for variant, result in results.items():
        solution = trajectory(result, initial[variant], times=times, rate=rate)
        color, style = ('#087e8b', '-') if variant == 'guarded' else ('#b7463a', '--')
        for ax, name in zip(axes, ('Z', 'variance', 'blocked')):
            ax.plot(times, solution[name], color=color, ls=style, lw=2,
                    label='With NAC' if variant == 'guarded' else 'Without NAC')
            data[variant+'_'+name] = solution[name]
    for ax, title in zip(axes, ('Expected chlorine count', 'Chlorine-count variance', 'Blocked with free sites')):
        ax.set(xlabel='Time', title=title)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=8)
    axes[0].set_ylim(-.05, 4.1)
    axes[1].set_ylim(bottom=-.02)
    axes[2].set_ylim(-.025, 1.025)
    fig.suptitle(f'Generated event model · κ={rate:g} · initial ring mask={mask}', fontsize=11)
    return fig, data


def context_nac_sliders(results, baseline, *, mask=0, rate=1., end=3.):
    import ipywidgets as w
    controls = {
        'rate': w.FloatSlider(description='κ', value=rate, min=0., max=max(2., rate),
                            step=.1, continuous_update=False, layout=w.Layout(width='45%')),
        'end': w.FloatSlider(description='End time', value=end, min=.5, max=max(5., end),
                           step=.5, continuous_update=False, layout=w.Layout(width='45%'))}
    return _live(controls, lambda **values: context_nac_curve(
        results, mask=mask, initial=baseline['initial'], **values), baseline)
