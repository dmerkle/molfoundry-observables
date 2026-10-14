#!/usr/bin/env python3
"""Time a one-carbon identity query using only public MolFoundry APIs.

Run: python scripts/benchmark_molfoundry_counts.py --carbons 8 24 54
Each measurement uses a fresh process. No mfobs imports, caching or JS access.
"""
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter


def measure(carbons, method):
    import molfoundry as mf

    # An ordinary linear alkane, C_n H_(2n+2), with explicit hydrogens.
    nodes = {i: 'C' for i in range(carbons)}
    edges = [(i, i+1) for i in range(carbons-1)]
    for carbon in range(carbons):
        for _ in range(4 - int(carbon > 0) - int(carbon < carbons-1)):
            hydrogen = len(nodes)
            nodes[hydrogen] = 'H'
            edges.append((carbon, hydrogen))
    gml = 'graph [\n' + '\n'.join(
        [f'node [ id {i} label "{label}" ]' for i, label in nodes.items()] +
        [f'edge [ source {u} target {v} label "-" ]' for u, v in edges]) + '\n]'
    record = {'molfoundry': version('molfoundry'), 'method': method, 'carbons': carbons,
              'host_vertices': len(nodes), 'expected_matches': carbons}
    start = perf_counter()
    host = mf.Graph.fromGMLString(gml, add=False)
    record['graph_seconds'] = perf_counter() - start
    start = perf_counter()
    rule = mf.Rule.fromGMLString(
        'rule [ context [ node [ id 0 label "C" ] ] ]', add=False)
    if method == 'direct':
        record['setup_seconds'] = perf_counter() - start
        print(json.dumps({**record, 'phase': 'before_direct_count'}), flush=True)
        start = perf_counter()
        record['matches'] = rule.numMatches([host], mode=mf.ApplyMode.UseAny)
        record['numMatches_seconds'] = perf_counter() - start
        assert record['matches'] == carbons, record
        print(json.dumps({**record, 'phase': 'complete'}), flush=True)
        return
    dg = mf.DG(graphDatabase=[host])
    record['setup_seconds'] = perf_counter() - start
    print(json.dumps({**record, 'phase': 'before_apply'}), flush=True)
    start = perf_counter()
    with dg.build() as builder:
        derivations = builder.apply([host], rule, mode=mf.ApplyMode.UseAny)
    record['apply_seconds'] = perf_counter() - start
    print(json.dumps({**record, 'phase': 'before_numMatches'}), flush=True)
    start = perf_counter()
    record['matches'] = sum(edge.numMatches(rule) for edge in derivations
                            if edge.numSources == 1)
    record['numMatches_seconds'] = perf_counter() - start
    assert record['matches'] == carbons, record
    print(json.dumps({**record, 'phase': 'complete'}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--carbons', nargs='+', type=int, default=[8, 24, 54])
    parser.add_argument('--timeout', type=float, default=30,
                        help='Maximum seconds per fresh process, including import')
    parser.add_argument('--worker', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--method', choices=['direct', 'dg'], default='direct',
                        help='Direct Rule.numMatches or DG.apply plus hyperedge counts')
    args = parser.parse_args()
    if any(n < 1 for n in args.carbons) or args.timeout <= 0:
        parser.error('carbon counts and timeout must be positive')
    if args.worker is not None:
        measure(args.worker, args.method)
        return
    for carbons in args.carbons:
        command = [sys.executable, str(Path(__file__).resolve()), '--worker', str(carbons),
                   '--method', args.method]
        try:
            run = subprocess.run(command, capture_output=True, text=True,
                                 timeout=args.timeout, check=True)
            print(run.stdout.splitlines()[-1], flush=True)
        except subprocess.TimeoutExpired as exc:
            lines = (exc.stdout or b'').decode().splitlines()
            partial = json.loads(lines[-1]) if lines else {'carbons': carbons}
            print(json.dumps({**partial, 'timeout_seconds': args.timeout}), flush=True)


if __name__ == '__main__':
    main()
