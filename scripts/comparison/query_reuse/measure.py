#!/usr/bin/env python3
"""Repeat sole-query discovery with one unchanged path backend and source rules.

A small reproducible query-reuse experiment, not a comparison with another
method. Each observation runs in a fresh serial process. The standard analysis is retained in raw JSON. Independent reference checks
run separately after the timed call.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT
sys.path.insert(0, str(PACKAGE / 'src'))
from mfobs.core.source_identity import require_source_identity

MOTIFS = ('AB', 'ABA', 'ABBA', 'AAAA', 'ABABA')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def worker(motif, destination):
    sys.path.insert(0, str(PACKAGE/'src'))
    from mfobs import ABA
    import molfoundry
    started = time.perf_counter()
    backend = ABA.motif(motif, gated=False)
    inputs = backend.specification()
    result = backend.jump_closure(inputs, linear_limit=24,
                     support_pattern_limit=64, support_expansion_limit=64,
                     diagnostic_timings=True)
    elapsed = time.perf_counter()-started
    result['linear'] = backend.linear_closure(inputs, maximum_dimension=24)
    result['linear']['input_definition_sha256'] = result['contract']['definition_sha256']
    # Darwin reports bytes, Linux KiB. This is whole-process high-water RSS,
    # including imported runtime; it is not the analysis's incremental memory.
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = rss if sys.platform == 'darwin' else rss*1024
    destination.write_text(json.dumps(result, indent=2)+'\n', encoding="utf-8")
    support = result['support_discovery']
    phases = result['diagnostics']['seconds']
    from mfobs.checks.terminal_paths import validate_result
    metadata = {
        'motif': motif, 'support_status': support['status'],
        'coordinates': list(support['patterns']),
        'suffixes': sorted(n for n,v in support['patterns'].items() if v['kind']=='suffix'),
        'event_status': result['pjc']['status'], 'events': result['pjc']['event_class_count'],
        'linear_status': result['linear']['status'], 'linear_dimension':len(result['linear']['basis']),
        'oracle':validate_result(result),
        'raw_composition_terms':sum(order['public_raw'] for row in support['trace']
                                   for order in row.get('composition_orders',{}).values()),
        'source_rules_sha256': hashlib.sha256(json.dumps(support['source_rules'],sort_keys=True).encode()).hexdigest(),
        'analysis_call_seconds': elapsed, 'phase_seconds': phases,
        'discovery_and_events_seconds': phases['support_discovery']+phases['event_compilation'],
        'peak_process_rss_bytes':rss_bytes, 'raw_sha256':sha(destination),
        'kernel_version':importlib.metadata.version('molfoundry'),
        'kernel_files':{p.relative_to(Path(molfoundry.__file__).resolve().parent).as_posix():sha(p)
                        for p in sorted(Path(molfoundry.__file__).resolve().parent.rglob('*'))
                        if p.is_file() and p.suffix in ('.py','.js','.so','.wasm')},
    }
    destination.with_suffix('.measurement.json').write_text(json.dumps(metadata,indent=2)+'\n', encoding="utf-8")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--worker', choices=MOTIFS)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--cpu-description', default='not recorded')
    args=parser.parse_args()
    args.out=args.out.resolve()
    if args.worker:
        worker(args.worker,args.out)
        return
    if args.repetitions<1: parser.error('positive repetitions required')
    identity = require_source_identity(ROOT)
    args.out.mkdir(parents=True,exist_ok=False)
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in sorted((PACKAGE/'src').rglob('*.py'))}
    records=[]
    for repetition in range(args.repetitions):
        # Deterministically rotate order between repetitions.
        order=MOTIFS[repetition%len(MOTIFS):]+MOTIFS[:repetition%len(MOTIFS)]
        for motif in order:
            destination=args.out/f'{motif}-{repetition+1}.json'
            command=[sys.executable,str(Path(__file__).resolve()),'--worker',motif,'--out',str(destination)]
            start=time.perf_counter()
            with destination.with_suffix('.log').open('w', encoding="utf-8") as log:
                try:
                    result=subprocess.run(command,cwd=PACKAGE,stdout=log,stderr=subprocess.STDOUT,timeout=60)
                except subprocess.TimeoutExpired:
                    raise RuntimeError(f'60-second wall limit exceeded: {motif}, repetition {repetition+1}')
            if result.returncode: raise RuntimeError(f'worker failed: {destination.with_suffix(".log")}')
            row=json.loads(destination.with_suffix('.measurement.json').read_text(encoding="utf-8"))
            row.update(repetition=repetition+1,fresh_process_wall_seconds=time.perf_counter()-start)
            records.append(row)
            print(motif,repetition+1,row['coordinates'],row['events'],round(row['analysis_call_seconds'],3),flush=True)
    assert sources=={p.relative_to(ROOT).as_posix():sha(p) for p in sorted((PACKAGE/'src').rglob('*.py'))}
    assert all(row['kernel_files']==records[0]['kernel_files'] for row in records)
    assert len({row['source_rules_sha256'] for row in records})==1, 'source rules changed with the query'
    summary={}
    for motif in MOTIFS:
        rows=[row for row in records if row['motif']==motif]
        for key in ('coordinates','suffixes','raw_composition_terms','events','linear_status','linear_dimension'):
            assert all(row[key]==rows[0][key] for row in rows), (motif,key)
        assert all(row['support_status']=='POLYNOMIAL_SUPPORT_CLOSED' and row['event_status']=='CLOSED'
                   and row['oracle']['status']=='PASS' for row in rows)
        summary[motif]={key:rows[0][key] for key in ('coordinates','suffixes','raw_composition_terms','events','linear_status','linear_dimension','oracle')}
        for key in ('analysis_call_seconds','discovery_and_events_seconds','peak_process_rss_bytes','fresh_process_wall_seconds'):
            summary[motif][key]={'median':statistics.median(row[key] for row in rows),
                                 'min':min(row[key] for row in rows),'max':max(row[key] for row in rows)}
    assert require_source_identity(ROOT) == identity, 'source identity changed during measurements'
    manifest={'status':'PASS','source_commit':identity['commit'], 'source_identity':identity,
              'script_sha256':sha(Path(__file__)), 'sources':sources,
              'platform':platform.platform(),'machine':platform.machine(),'python':sys.version,
              'host':platform.node(),
              'cpu':args.cpu_description,
              'kernel_version':records[0]['kernel_version'],'kernel_files':records[0]['kernel_files'],
              'protocol':{'repetitions':args.repetitions,'serial':True,'fresh_process':True,'wall_limit_seconds':60,
                          'linear_limit':24,'support_pattern_limit':64,'support_expansion_limit':64,'gated':False,
                          'timing_scope':'instrumented Algorithm 2 discovery+events; analysis call also includes EMGF export, but excludes Algorithm 1, reference checks, imports and JSON writing',
                          'memory_scope':'whole-process high-water RSS including imported runtime',
                          'term_scope':'sum of returned raw composition products in both orders, before family filtering and signed cancellation; not candidates visited'},
              'source_rules_sha256':records[0]['source_rules_sha256'],'cases':summary,'observations':records}
    (args.out/'summary.json').write_text(json.dumps(manifest,indent=2)+'\n', encoding="utf-8")
    print('PASS:',args.out/'summary.json')

if __name__=='__main__':main()
