#!/usr/bin/env python3
"""Recompute four chemical stages and one explicit budget outcome in fresh processes."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=int, default=300)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cases = [
        ('baseline', ['--preset', 'baseline']),
        ('noedge', ['--preset', 'noedge']),
        ('aba', ['--motif', 'ABA']),
        ('aba-gated', ['--motif', 'ABA', '--gate']),
        ('limited-aba', ['--motif', 'ABA', '--molecular-overlap-limit', '0']),
    ]
    environment = dict(os.environ, PYTHONPATH=str(root / 'src'))
    summary = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        'python': platform.python_version(), 'platform': platform.platform(),
        'wall_clock_limit_per_case_seconds': args.seconds,
        'scope': 'fresh molecular route only; neither frozen artifact nor historical typed benchmarks rerun',
        'cases': {},
    }
    for name, options in cases:
        path = output / (name + '.json')
        command = [sys.executable, 'run.py', 'rop-analysis', *options,
                   '--linear-limit', '24',
                   '--diagnostic-timings', '--out', str(path)]
        log = output / (name + '.log')
        with log.open('w', encoding='utf-8') as stream:
            try:
                subprocess.run(command, cwd=root, env=environment, stdout=stream,
                               stderr=subprocess.STDOUT, timeout=args.seconds, check=True)
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
                summary['cases'][name] = {'status': 'FAILED_RUN', 'reason': str(exc), 'command': command}
                (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
                raise
        data = json.loads(path.read_text(encoding='utf-8'))
        from mfobs.checks.molecular_rop import validate_result
        checks = validate_result(data)
        expected = 'UNKNOWN_COMPOSITION_LIMIT' if name == 'limited-aba' else 'ANALYZED'
        if data['status'] != expected:
            raise AssertionError((name, data['status'], expected))
        folder = output / name
        folder.mkdir()
        for identifier, source in data['support_discovery']['source_rules'].items():
            (folder / ('rule-' + identifier + '.gml')).write_text(source['gml'], encoding='utf-8')
        for identifier, pattern in data['support_discovery']['patterns'].items():
            (folder / ('observable-' + identifier + '.gml')).write_text(pattern['gml'], encoding='utf-8')
        summary['cases'][name] = {
            'command': command, 'status': data['status'],
            'selected_queries': data['support_discovery']['selected_queries'],
            'coordinates': list(data['support_discovery']['patterns']),
            'support': data['support_discovery']['status'],
            'events': data['support_discovery']['events']['status'],
            'event_count': data['pjc'].get('event_class_count'),
            'linear': data['linear']['status'], 'linear_dimension': len(data['linear'].get('basis', [])),
            'emgf_order': data['pjc'].get('operator', {}).get('maximum_differential_order'),
            'molecular_validation': checks['molecular_graph_validation'],
            'typed_comparison': checks['typed_reference_comparison'],
            'diagnostics': data['diagnostics'],
            'analysis_sha256': sha256(path.read_bytes()).hexdigest(),
        }
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        print(name, data['status'], data['linear']['status'], data['pjc'].get('event_class_count'), flush=True)
    manifest = {p.relative_to(output).as_posix(): sha256(p.read_bytes()).hexdigest()
                for p in sorted(output.rglob('*')) if p.is_file()}
    (output / 'SHA256.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
