"""Command-line entry points. Commands read actual inputs; no canned ODE output."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction
import argparse,csv,hashlib,importlib.metadata,json,platform,sys,traceback
from mfobs.inputs import load,InputError


def parser():
    p=argparse.ArgumentParser(prog='mf-observables',description='Conditional-DPO observable dynamics with exact invariant counting.')
    p.add_argument('--debug',action='store_true',help='Print a Python traceback on failure')
    s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('inspect',help='Show exact input files, rates, assumptions and initial observables');q.add_argument('model')
    q=s.add_parser('run',help='Recompute compositions, reduce, and print/export the ODE');q.add_argument('model')
    q.add_argument('--out',default=None);q.add_argument('--mode',choices=['restricted','public'])
    q.add_argument('--progress',action='store_true');q.add_argument('--max-expansions',type=int,default=1000)
    q.add_argument('--force',action='store_true',help='Replace known output files from a previous run')
    q=s.add_parser('ode',help='Print a previously computed ODE, without recomposition');q.add_argument('result');q.add_argument('--normalized',action='store_true')
    q=s.add_parser('compress',help='Exact seed-preserving rate-family linear compression of a closed result')
    q.add_argument('result');q.add_argument('--out',required=True)
    q.add_argument('--force',action='store_true',help='Replace known output files, never the source result')
    q=s.add_parser('verify',help='Check computed coefficients and molecular/finite-state semantics');q.add_argument('model');q.add_argument('--result',required=True)
    q.add_argument('--full',action='store_true',help='Use full molecular audit rather than the small smoke grid')
    q.add_argument('--multiplicities',action='store_true',help='Compare both composition backends on all production pairs')
    q.add_argument('--conditions',action='store_true',help='Independent adjacency-condition transport diagnostics')
    q.add_argument('--out',default=None)
    q=s.add_parser('solve',help='Numerically evaluate the EXPORTED linear ODE; requires scipy');q.add_argument('result')
    q.add_argument('--parameters',required=True);q.add_argument('--initial');q.add_argument('--until',type=float,default=10)
    q.add_argument('--points',type=int,default=101);q.add_argument('--out',required=True)
    q=s.add_parser('derive',help='Generic single observable/rule commutator; no example-specific reducer')
    q.add_argument('--observable',required=True);q.add_argument('--rule',required=True);q.add_argument('--out',required=True)
    q.add_argument('--connected',action='store_true');q.add_argument('--include-empty',action='store_true')
    q=s.add_parser('rop-analysis',help='Analyze one preset or a custom A/B motif query without rebuilding reports')
    q.add_argument('--preset',choices=['baseline','noedge','copolymer','combined'])
    q.add_argument('--alphabet',default=None,help='Comma-separated subset of the literal chemical symbols A,B (custom mode only)')
    q.add_argument('--motif',action='append',default=[],help='Selected sequence motif, repeatable; for example --motif ABBA')
    q.add_argument('--gate',action='store_true',help='Enable the existing active NoEdge(Pd,N) gate (custom mode only)')
    q.add_argument('--linear-limit',type=int,default=24,help='Maximum retained dimension for linear discovery')
    q.add_argument('--molecular-overlap-limit',type=int,default=2000000,help='Per-order candidate overlap budget for the direct molecular route')
    q.add_argument('--support-pattern-limit',type=int,default=64,help='Maximum retained patterns for support discovery')
    q.add_argument('--support-expansion-limit',type=int,default=64,help='Maximum expanded patterns for support discovery; zero returns UNKNOWN')
    q.add_argument('--diagnostic-timings',action='store_true',help='Record wall-clock analysis phase timings in JSON; does not impose a timeout')
    q.add_argument('--out',required=True,help='JSON analysis output; this command never writes a public report')
    s.add_parser('doctor',help='Report runtime, dependency constraints and engine fingerprint')
    return p


def environment():
    versions={}
    for name in ['molfoundry','sympy','networkx','requests','numpy','scipy','rdkit','pytest']:
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]=None
    warnings=[]
    try:
        from packaging.version import Version
        if versions['requests'] and Version(versions['requests'])<Version('2.33'):
            warnings.append('Bundled MolFoundry declares requests>=2.33; this environment has an older Requests. Local computations can be tested, but this is not a dependency-clean or network-functionality validation.')
    except ImportError:warnings.append('packaging not installed; version constraint not checked')
    binaries={}
    try:
        importlib.import_module('molfoundry')
        dist=importlib.metadata.distribution('molfoundry')
        for rel in dist.files or []:
            f=Path(dist.locate_file(rel))
            if f.is_file() and (f.suffix in ['.so','.dylib','.pyd','.js'] or str(rel)=='molfoundry/_engine.py'):
                binaries[str(rel)]=hashlib.sha256(f.read_bytes()).hexdigest()
        engine_loaded=True
    except Exception as exc:engine_loaded=False;warnings.append(f'Engine import failed: {type(exc).__name__}: {exc}')
    return {'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),
            'versions':versions,'engine_loaded':engine_loaded,'binary_sha256':binaries,'warnings':warnings}


def inspect(inputs):
    import molfoundry as mf

    d=inputs.data
    print(f"{d['title']}\nInput: {inputs.path}\nAdapter: {d['adapter']}\n")
    print('Rules (rate PER INJECTIVE MATCH):')
    for r in d['rules']:
        rule=mf.Rule.fromGMLString(inputs.file(r['file']).read_text(encoding='utf-8'), add=False)
        print(f"  {r['id']}: {r['per_match']} * {r['rate']} [{r['kind']}]  {r['file']}")
        print(f"    conditions: NoEdge={len(rule.matchConstraintsNoEdge)} "
              f"(created-edge={sum(c.implied for c in rule.matchConstraintsNoEdge)}), "
              f"adjacency={len(rule.matchConstraintsAdjacency)}, "
              f"nested={'yes' if rule.matchCondition is not None else 'no'}")
    print('\nInitially queued: '+', '.join(d['initial_observables']))
    print('Available input pattern files: '+', '.join(o['id'] for o in d['observables']))
    print('Composition: '+json.dumps(d['composition']))
    print('Forbidden family: '+str(d['forbidden']))
    if d.get('specifications'):
        print('Specification layers: '+json.dumps(d['specifications']))
        print('Chemical refinement: '+json.dumps(inputs.specification_certificate()))
    print('State language: '+d['state_language']);print(json.dumps(inputs.read('state_language'),indent=2,ensure_ascii=False))
    print('\nReduction input: '+d['reduction']);print(json.dumps(inputs.read('reduction'),indent=2,ensure_ascii=False))
    print('\nInitial state (only for exported initial values): '+json.dumps(d['initial_state']))
    if d.get('initial_mixture'):
        print('Concrete initial component mixture: '+d['initial_mixture'])
        print(json.dumps(inputs.initial_mixture(),indent=2,ensure_ascii=False))
    print('Numerical rates (NOT used in symbolic synthesis): '+d['parameters'])
    print('Reference result is NEVER read by run; only by verify.')
    print('Definition SHA256: '+inputs.definition_digest())


def run_command(a):
    inputs=load(a.model)
    if a.mode:
        inputs.data['composition']['mode']=a.mode
        if a.mode.startswith('restricted') and not inputs.data['forbidden']:raise InputError('This model does not supply a forbidden family')
    out=Path(a.out or ('out/'+inputs.data['id']))
    if (out/'result.json').exists() and not a.force:raise InputError(f'{out}/result.json already exists; use --force to replace this output')
    out.mkdir(parents=True,exist_ok=True)
    from mfobs.core.worklist import run
    from mfobs.outputs import export,ode_text
    result=run(inputs,progress=a.progress,max_expansions=a.max_expansions,checkpoint=out/'checkpoint.json')
    result['runtime']=environment()
    export(result,out)
    (out/'inputs.json').write_text(json.dumps(inputs.data,indent=2)+'\n', encoding='utf-8')
    print(ode_text(result))
    print('Rounds: '+' -> '.join(str(r['after']) for r in result['rounds']))
    print(f"Raw compositions: {result['raw_compositions']}; ordered calls: {result['ordered_composition_calls']}")
    files='ode.txt, ode.json, model.sbml.xml, ode.tex, result.json, trace.json, matrix.csv, observables/, identity-patterns/' if result['closed'] else 'result.json, trace.json, checkpoint.json, observables/, identity-patterns/; no closed ODE matrix'
    print(f"Outputs: {out.resolve()} ({files})")
    return 0 if result['closed'] else 2


def verify_command(a):
    inputs=load(a.model);result=json.loads(Path(a.result).read_text(encoding='utf-8'))
    from mfobs.checks.reference import verify
    report={'model':inputs.data['id'],'reference':verify(inputs,result)}
    if inputs.data['adapter']=='mfobs.models.chain_mixture:Adapter':
        from mfobs.checks.invariant_discovery import verify as invariant_verify
        report['graphs']=invariant_verify(inputs,result,full=a.full)
    elif inputs.data['id']=='ring8':
        from mfobs.checks.ring import verify as ring_verify
        report['graphs']=ring_verify(inputs,result)
    else:report['graphs']={'performed':False,'reason':'no independent model oracle installed'}
    if a.multiplicities:
        from mfobs.checks.multiplicities import verify as counter_verify
        report['multiplicities']=counter_verify(inputs,result,full=True)
    if a.conditions:
        from mfobs.checks.conditions import verify as conditions_verify
        report['condition_transport']=conditions_verify()
    report['passed']=True;report['runtime']=environment()
    out=Path(a.out) if a.out else Path(a.result).parent/'verification.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n', encoding='utf-8')
    print(f'PASS: {out.resolve()}')
    return 0


def compress_command(a):
    from mfobs.invariants.compress import compress
    from mfobs.outputs import export
    source=Path(a.result).resolve();out=Path(a.out).resolve()
    if (out/'result.json').resolve()==source:
        raise InputError('Compression requires a separate output directory; never overwrite its source result')
    if (out/'result.json').exists() and not a.force:
        raise InputError(f'{out}/result.json already exists; use --force to replace this output')
    result=compress(json.loads(source.read_text(encoding='utf-8')))
    result['runtime']=environment()
    export(result,out)
    print(f"Exact rate-family compression: {result['source_dimension']} -> {result['dimension']} coordinates.")
    print('Coordinates are rational linear combinations of source observables, not individual graph patterns.')
    print(f'Outputs: {out} (ODE, transformation.csv, basis-combinations.json, seed-reconstruction.json, result.json)')
    return 0


def solve_command(a):
    import numpy as np
    import sympy as sp
    from scipy.linalg import expm
    from mfobs.outputs import matrix
    d=json.loads(Path(a.result).read_text(encoding='utf-8'))
    if not d['closed']:raise InputError('Cannot solve an UNKNOWN/partial output')
    params=json.loads(Path(a.parameters).read_text(encoding='utf-8'));required={r['rate'] for r in d['rates']}
    if set(params)!=required:raise InputError(f'Parameters must be exactly {sorted(required)}')
    values={k:Fraction(v) for k,v in params.items()}
    if any(v<0 for v in values.values()):raise InputError('Rates must be nonnegative')
    initial=json.loads(Path(a.initial).read_text(encoding='utf-8')) if a.initial else d['initial_values']
    if set(initial)!=set(d['basis']):raise InputError('Initial expectation vector must cover the entire output basis in its exported coordinates')
    if not 2<=a.points<=100000 or not 0<=a.until<float('inf'):raise InputError('points must be 2..100000 and until finite/nonnegative')
    am=matrix(d).subs({sp.Symbol(k):sp.Rational(v.numerator,v.denominator) for k,v in values.items()})
    aa=np.array(am.tolist(),dtype=float);y0=np.array([float(Fraction(initial[k])) for k in d['basis']])
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='', encoding='utf-8') as f:
        wr=csv.writer(f);wr.writerow(['t',*['mu_'+k for k in d['basis']]])
        for t in np.linspace(0,a.until,a.points):wr.writerow([float(t),*(expm(float(t)*aa)@y0).tolist()])
    print(f'Wrote {out.resolve()}; numerical evaluation of the exported exact model ODE, not a simulation or experimental fit.')
    return 0


def derive_command(a):
    from mfobs.core.engine import mf,commutator
    obs=mf.Rule.fromGMLString(Path(a.observable).read_text(encoding='utf-8'),add=False)
    rule=mf.Rule.fromGMLString(Path(a.rule).read_text(encoding='utf-8'),add=False)
    c=commutator(obs,rule,connected=a.connected,include_empty=a.include_empty)
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True);terms=[]
    for i,(key,co) in enumerate(sorted(c.observable_coefficients.items())):
        filename=f'jump-{i:04}.gml';(out/filename).write_text(c.observable_representatives[key].getGMLString(), encoding='utf-8')
        terms.append({'canonical_key':key,'coefficient':co,'file':filename})
    report={'scope':'single commutator, no model-specific reduction or closure claim','left_raw':c.left_raw,'right_raw':c.right_raw,'terms':terms}
    (out/'commutator.json').write_text(json.dumps(report,indent=2)+'\n', encoding='utf-8')
    print(f'{len(terms)} jump terms, {c.left_raw+c.right_raw} raw compositions. See {out}/commutator.json')
    return 0


def rop_analysis_command(a):
    """Run the existing composition pipeline for a configurable query only."""

    from mfobs.api import analyze as analyze_model
    from mfobs.applications import MolecularROP
    from mfobs.models.molecular_rop import inputs as rop

    custom_options = bool(a.motif) or a.alphabet is not None or a.gate
    if a.preset:
        if custom_options:
            raise InputError(
                '--preset cannot be combined with --motif, --alphabet, or --gate; '
                'omit --preset for a custom analysis'
            )
        selection = a.preset
    else:
        if not a.motif:
            raise InputError(
                'custom ROP analysis needs at least one --motif, for example --motif ABBA'
            )
        if a.alphabet is None:
            raw_alphabet = 'A,B'
        else:
            raw_alphabet = a.alphabet
        alphabet = tuple(
            value.strip() for value in raw_alphabet.split(',')
            if value.strip()
        )
        if not alphabet:
            raise InputError('--alphabet was provided but names no chemical symbols; use A, B, or A,B')
        selection = rop.configured_variant(
            motifs=tuple(a.motif), gated=a.gate, alphabet=alphabet,
        )
    backend = MolecularROP(selection)
    analysis = analyze_model(
        backend.specification(), backend=backend,
        linear_limit=a.linear_limit,
        diagnostic_timings=a.diagnostic_timings,
        support_pattern_limit=a.support_pattern_limit,
        support_expansion_limit=a.support_expansion_limit,
        molecular_overlap_limit=a.molecular_overlap_limit,
    ).record()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(analysis, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    configuration = analysis['configuration']
    if analysis.get('status', '').startswith('UNKNOWN'):
        print(
            f"{analysis['status']}: {out.resolve()} · "
            f"support={analysis['support_discovery']['status']} · "
            f"events={analysis['support_discovery']['events']['status']} · "
            "incomplete obligations and trace retained in JSON"
        )
        return 0
    print(
        f"PASS: {out.resolve()} · {configuration['model']['identifier']} · "
        f"events={len(analysis['composition']['event_classes'])} · "
        f"linear={analysis['linear']['status']} · pJC={analysis['pjc']['status']} · "
        f"configuration {configuration['sha256']}"
    )
    print(
        "Scope: exact analysis is relative to the declared ROP state language "
        "and its explicit model assumptions; details are in the JSON record."
    )
    return 0


def main(argv=None):
    a=parser().parse_args(argv)
    try:
        if a.command=='doctor':
            d=environment();print(json.dumps(d,indent=2));return 0 if d['engine_loaded'] else 1
        if a.command=='inspect':inspect(load(a.model));return 0
        if a.command=='run':return run_command(a)
        if a.command=='ode':
            from mfobs.outputs import ode_text
            print(ode_text(json.loads(Path(a.result).read_text(encoding='utf-8')),a.normalized),end='');return 0
        if a.command=='verify':return verify_command(a)
        if a.command=='compress':return compress_command(a)
        if a.command=='solve':return solve_command(a)
        if a.command=='derive':return derive_command(a)
        if a.command=='rop-analysis':return rop_analysis_command(a)
    except (Exception,KeyboardInterrupt) as exc:
        if a.debug:traceback.print_exc()
        else:print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
