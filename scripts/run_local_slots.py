#!/usr/bin/env python3
"""Bounded non-enumerative local-slot synthesis. Old Ring8 outputs are untouched."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mfobs.inputs import load
from mfobs.core.worklist import run
from mfobs.outputs import export


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',choices=['local','degree2','no-coop','all'],default='all')
    parser.add_argument('--max-expansions',type=int,default=6)
    args=parser.parse_args()
    variants=['no-coop','local','degree2'] if args.model=='all' else [args.model]
    for tag in variants:
        inputs=load(ROOT/f'examples/local-slots/model-{tag}.json')
        result=run(inputs,progress=True,max_expansions=args.max_expansions)
        result['source_sha256']={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in [
            Path(__file__),ROOT/'src/mfobs/invariants/local_slots.py',ROOT/'src/mfobs/models/local_slots.py']}
        out=ROOT/f'out/local-slots/{tag}'
        export(result,out)
        print(json.dumps({'model':tag,'closed':result['closed'],'patterns':len(result['basis']),
            'expanded':len(result['steps']),'raw_compositions':result['raw_compositions'],
            'local_reductions':result['local_reduction_statistics'],'seconds':result['seconds'],
            'output':str(out)},indent=2),flush=True)


if __name__=='__main__':main()
