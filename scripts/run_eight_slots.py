#!/usr/bin/env python3
"""Discover the all-topology eight-slot ODE without a finite-state oracle."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from mfobs import EightSlots, analyze


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-expansions', type=int, default=100)
    parser.add_argument('--out', type=Path, default=ROOT/'out/eight-slots')
    args = parser.parse_args()
    out = args.out.resolve()
    backend = EightSlots()
    result = analyze(backend.specification(), backend=backend, output=out, max_expansions=args.max_expansions, progress=True).record()
    print(json.dumps({'status': result['status'], 'patterns': len(result['basis']),
        'expanded': len(result['steps']), 'pending': result['pending'],
        'unknown': result['unknown'], 'raw_compositions': result['raw_compositions'],
        'seconds': result['seconds'], 'output': str(out)}, indent=2), flush=True)


if __name__ == '__main__':
    main()
