#!/usr/bin/env python3
"""Execute numbered notebooks in fresh kernels and bind saved outputs to their cells."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
CHAPTERS = ROOT / 'notebooks/examples'
sys.path.insert(0, str(ROOT / 'notebooks'))
from execute import prepare_saved_outputs
from mfobs.provenance import runtime_fingerprint
from mfobs.core.source_identity import source_identity


def source_digest(book):
    return sha256(json.dumps([(c.cell_type, c.source) for c in book.cells],
                            ensure_ascii=False).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('chapters', nargs='*', help='Stems such as 02_aba; default: all numbered chapters')
    args = parser.parse_args()
    paths = [CHAPTERS/(name+'.ipynb') for name in args.chapters] if args.chapters else sorted(CHAPTERS.glob('[0-9][0-9]_*.ipynb'))
    if not args.chapters:
        paths += sorted((ROOT / 'notebooks/visual/comparisons').glob('*.ipynb'))
    revision = source_identity(ROOT)['commit']
    for path in paths:
        book = nbformat.read(path, as_version=4)
        digest = source_digest(book)
        print('Executing', path.name, flush=True)
        started = time.monotonic()
        NotebookClient(book, timeout=2400, kernel_name='molfoundry',
                       resources={'metadata': {'path': str(path.parent)}}).execute()
        if source_digest(nbformat.read(path, as_version=4)) != digest:
            raise RuntimeError(f'{path}: edited during execution; refusing to overwrite')
        book.metadata['execution'] = {
            'runtime': runtime_fingerprint(),
            'cell_source_sha256': digest,
            'source_commit': revision,
            'method': 'nbclient; fresh molfoundry kernel',
        }
        # Live widget state is recreated by execution; static figure outputs remain.
        prepare_saved_outputs(book)
        nbformat.write(book, path)
        print(f'Passed {path.name} in {time.monotonic()-started:.1f}s', flush=True)


if __name__ == '__main__':
    main()
