#!/usr/bin/env python3
"""Execute the full notebook in a fresh kernel; retain partial output on error."""
from __future__ import annotations
import argparse
from pathlib import Path
import nbformat
from nbclient import NotebookClient


def prepare_saved_outputs(book):
    """Keep static figures, but never publish controls tied to a dead kernel."""
    book.metadata.pop('widgets', None)
    for cell in book.cells:
        for output in cell.get('outputs', []):
            data = output.get('data', {})
            if data.pop('application/vnd.jupyter.widget-view+json', None) is not None:
                data.setdefault('text/html',
                    '<p><em>Restart Kernel and Run All to activate the controls.</em></p>')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--notebook', required=True, type=Path, help='One example notebook (the index has no executable analysis)')
    parser.add_argument('--kernel', default='molfoundry')
    parser.add_argument('--timeout', type=int, default=900, help='Seconds per cell')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; select a new path')
    if args.timeout < 1:
        parser.error('timeout must be positive')
    source = args.notebook.resolve()
    if not source.is_file():
        parser.error(f'Notebook does not exist: {source}')
    nb = nbformat.read(source, as_version=4)
    if not any(cell.cell_type == 'code' for cell in nb.cells):
        parser.error('Choose an example notebook; this notebook is only an index')
    for cell in nb.cells:
        if cell.cell_type == 'code':
            cell.outputs = []
            cell.execution_count = None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        NotebookClient(nb, kernel_name=args.kernel, timeout=args.timeout,
                       resources={'metadata': {'path': str(source.parent)}},
                       allow_errors=False).execute()
    finally:
        prepare_saved_outputs(nb)
        with args.output.open('x', encoding='utf-8') as stream:
            nbformat.write(nb, stream)
    print(f'PASS: {args.output}')


if __name__ == '__main__':
    main()
