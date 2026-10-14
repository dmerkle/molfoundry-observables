"""Presentation and numerical consumers for the remaining visual examples.

No synthesis equations are supplied here. Exact curves consume recorded
coefficients/events; approximation curves call the existing reference models.
"""
from fractions import Fraction
from html import escape
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.integrate import solve_ivp
from IPython.display import display, HTML, Math
import paper_view as view
import inspection as files


def gml_table(patterns, title='Discovered patterns'):
    rows=[]
    for name, record in patterns.items():
        gml=record if isinstance(record,str) else record['gml']
        meaning='' if isinstance(record,str) else record.get('meaning',record.get('origin',''))
        if isinstance(meaning, dict):
            meaning = meaning.get('meaning', meaning.get('kind', '')).replace('-', ' ')
        rows.append([escape(name),files.gml_link(gml),escape(str(meaning))])
    display(HTML('<details><summary>'+escape(title)+'</summary>'+files.html_table(['Count','Pattern','Meaning'],rows)+'</details>'))


def eight_matrix(result,rates):
    names=result['basis']; index={n:i for i,n in enumerate(names)}
    if not result['closed'] or result['pending'] or result['unknown']:raise ValueError('No closed symbolic ODE')
    channels={c for row in result['equations_by_reaction'].values() for c in row}
    if set(rates)!=channels or any(not np.isfinite(v) or v<0 for v in rates.values()):raise ValueError('Supply all nonnegative channel rates')
    A=np.zeros((len(names),len(names)))
    for n,terms in result['equations_by_reaction'].items():
        for channel,row in terms.items():
            for target,c in row.items(): A[index[n],index[target]]+=rates[channel]*float(Fraction(c))
    return names,A


def eight_initial(result, bits):
    # Evaluation happens only after synthesis and cannot influence discovery.
    from mfobs.models.eight_slots import TOPOLOGIES, polynomial_count
    if len(bits)!=8 or any(x not in (0,1) for x in bits):raise ValueError('Supply eight binary occupancies')
    return {tag:{name:float(polynomial_count((result['basis_polynomials'][name] if 'basis_polynomials' in result else result['symbolic_patterns'][name]['polynomial']),bits,lengths))
                 for name in result['basis']} for tag,lengths in TOPOLOGIES.items()}


def topology_picture(bits):
    from mfobs.models.eight_slots import TOPOLOGIES
    fig,axes=plt.subplots(1,3,figsize=(12,2.8),layout='constrained')
    for ax,(tag,lengths) in zip(axes,TOPOLOGIES.items()):
        offset=0
        for component,n in enumerate(lengths):
            angles=np.arange(n)*2*np.pi/n
            x=np.cos(angles)+2.7*component;y=np.sin(angles)
            ax.plot(np.r_[x,x[0]],np.r_[y,y[0]],color='#7b8995')
            ax.scatter(x,y,c=bits[offset:offset+n],cmap='Blues',vmin=0,vmax=1,s=110,edgecolors='#264057',zorder=3)
            offset+=n
        ax.set_title(tag);ax.set_aspect('equal');ax.axis('off')
    fig.suptitle('C-next core · blue = occupied own P–L bound edge; white = absent')
    display(fig);plt.close(fig)


def eight_curves(result,initial,rates,*,end):
    names,A=eight_matrix(result,rates);times=np.linspace(0,end,151)
    values={tag:expm_multiply(csr_matrix(A),np.array([row[n] for n in names]),start=0,stop=end,num=len(times)) for tag,row in initial.items()}
    fig,axes=plt.subplots(1,2,figsize=(12,3.8),layout='constrained');col=names.index('P0')
    for tag,y in values.items():
        axes[0].plot(times,y[:,col],label=tag)
        if tag!='8':axes[1].plot(times,y[:,col]-values['8'][:,col],label=tag+' minus 8')
    axes[0].set(ylabel='Expected occupied slots',ylim=(0,8),title='Same symbolic ODE, different initial vectors')
    axes[1].set(ylabel='Difference',title='Topology comparison')
    for ax in axes:ax.set_xlabel('Time');ax.legend()
    display(fig);plt.close(fig)
    return {'time':times,'values':values,'rates':rates,'options':{'end':end}}


def rate_controls(callback, rates, options):
    """Changing rates keeps the original time/seed/run settings, explicitly."""
    import ipywidgets as widgets
    controls={name:widgets.FloatSlider(value=value,min=0,max=max(2.,value*3),step=.05,description=name,continuous_update=False)
              for name,value in rates.items()}
    output=widgets.Output();state={'rates':dict(rates),'options':dict(options),'result':None}
    button=widgets.Button(description='Update curves')
    def update(_):
        with output:
            output.clear_output(wait=True)
            state['rates']={name:control.value for name,control in controls.items()}
            state['result']=callback(state['rates'],**state['options'])
    button.on_click(update)
    display(widgets.VBox([*controls.values(),button,output]))
    return state


def coefficient_equations(result,names=None):
    for name,channels in result['equations_by_reaction'].items():
        if names is not None and name not in names:continue
        expression=sum(sp.Symbol(channel)*sum(sp.Rational(c)*sp.Symbol(target) for target,c in row.items()) for channel,row in channels.items())
        display(Math(r'\frac{d}{dt}\mathbb{E}['+name+r'] = '+sp.latex(expression)))


def walker_matrix(record,rates):
    names=['E','C'];symbols={n:sp.Symbol(n) for n in names+list(rates)}
    equations=record['comparison']['physical_equations']
    A=np.array([[float(sp.diff(sp.sympify(equations[n],locals=symbols),symbols[m]).subs({symbols[k]:v for k,v in rates.items()})) for m in names] for n in names])
    if np.max(np.abs(A.sum(axis=0)))>1e-12:raise ValueError('Walker probability is not conserved')
    return names,A


def walker_curve(record,rates,*,end,initial):
    names,A=walker_matrix(record,rates);time=np.linspace(0,end,151)
    values=expm_multiply(csr_matrix(A),np.array([initial[n] for n in names]),start=0,stop=end,num=len(time))
    fig,ax=plt.subplots(figsize=(10,3.8),layout='constrained')
    for i,name in enumerate(names):ax.plot(time,values[:,i],label=name)
    ax.set(xlabel='Time',ylabel='Configuration probability',title='DNA walker: configuration probabilities');ax.legend();display(fig);plt.close(fig)
    return {'time':time,'values':values,'names':names,'rates':rates}


def bruss_curve(result,initial,rates,*,end,paths,seed,volume):
    from mfobs.models.brusselator_growth_chemistry import deterministic_rhs, event_rates
    times=np.linspace(0,end,401)
    solution=solve_ivp(lambda t,y:deterministic_rhs(t,y,rates),(0,end),initial,t_eval=times,method='DOP853',rtol=1e-9,atol=1e-11)
    if not solution.success:raise RuntimeError(solution.message)
    count_rates=event_rates(rates,volume)
    counts={'X':round(initial[0]*volume),'Y':round(initial[1]*volume),'M':int(initial[2]),'N':int(initial[3])}
    names=view.event_system(result).names
    seeds=np.random.SeedSequence(seed).spawn(paths)
    samples=np.stack([view.simulate_events(result,counts,count_rates,times,np.random.default_rng(s),max_events=200000)[0] for s in seeds])
    x=samples[:,:,names.index('X')]/volume
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for y in x[:8]:axes[0].step(times,y,where='post',alpha=.5)
    axes[0].plot(times,solution.y[0],color='black',label='Macroscopic approximation');axes[0].legend(fontsize=8)
    axes[1].plot(solution.y[0],solution.y[1]);axes[1].set(xlabel='x',ylabel='y',title='Approximate phase portrait')
    axes[2].plot(times,x.mean(0),label='SSA mean')
    error=1.96*x.std(0,ddof=1)/np.sqrt(paths)
    axes[2].fill_between(times,x.mean(0)-error,x.mean(0)+error,alpha=.2,label='±1.96 MC SE');axes[2].legend(fontsize=8)
    for ax in (axes[0],axes[2]):ax.set(xlabel='Time',ylabel='X / V')
    axes[0].set_title('Computed-event sample paths');axes[2].set_title('SSA estimate, not a closed mean ODE')
    display(fig);plt.close(fig)
    return {'time':times,'approximate':solution.y,'samples':samples,'rates':rates,'options':dict(end=end,paths=paths,seed=seed,volume=volume)}


def growing_curve(initial,rates,*,end):
    from mfobs.models.structural_nac_oscillator import one_site_rhs
    times=np.linspace(0,end,501);values={}
    fig,axes=plt.subplots(1,2,figsize=(12,3.8),layout='constrained')
    for ax,g in zip(axes,(0.,rates['g'])):
        solution=solve_ivp(lambda t,y:one_site_rhs(y,rates['r'],g),(0,end),initial,t_eval=times,method='DOP853',rtol=1e-10,atol=1e-12)
        if not solution.success:raise RuntimeError(solution.message)
        values[str(g)]=solution.y
        for i,y in enumerate(solution.y):ax.plot(times,y,label=f'q{i}')
        ax.set(xlabel='Time',ylabel='Approximate fraction',title=f'R1 approximation · g={g:g}',ylim=(0,1));ax.legend()
    display(fig);plt.close(fig)
    return {'time':times,'values':values,'rates':rates,'options':{'end':end}}
