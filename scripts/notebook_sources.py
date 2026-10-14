#!/usr/bin/env python3
"""Export percent-cell Python from the authoritative numbered notebooks."""
from pathlib import Path
import argparse
import json
import re
from html import unescape

ROOT = Path(__file__).resolve().parents[1]


def source(book):
    parts = ['# Generated from the adjacent .ipynb by scripts/notebook_sources.py.\n# Edit the notebook, then regenerate this file.\n']
    for cell in book['cells']:
        text = ''.join(cell['source']).rstrip()
        if cell['cell_type'] == 'markdown':
            text = re.sub(r'<nav\b.*?</nav>', '', text, flags=re.S)
            text = re.sub(r'<h2\b[^>]*>(.*?)</h2>', lambda m: '## ' + unescape(m[1]), text)
            text = re.sub(r'<a id="[^"]*"></a>', '', text).strip()
            text = re.sub(r'<p>(<strong.*?</strong>.*?)</p>',
                          lambda m: '## ' + re.sub(r'<[^>]+>', '', m[1]), text)
            text = unescape(text)
            if not text:
                continue
            parts.append('# %% [markdown]\n' + '\n'.join('# ' + line if line else '#' for line in text.splitlines()) + '\n')
        elif cell['cell_type'] == 'code':
            parts.append('# %%\n' + text + '\n')
    return '\n'.join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    paths = [ROOT/'notebooks/paper_examples.ipynb', *sorted((ROOT/'notebooks/examples').glob('[0-9][0-9]_*.ipynb')), *sorted((ROOT/'notebooks/visual/comparisons').glob('*.ipynb'))]
    for path in paths:
        expected = source(json.loads(path.read_text(encoding='utf-8')))
        target = path.with_suffix('.py')
        if args.check:
            if not target.exists() or target.read_text(encoding='utf-8') != expected:
                raise SystemExit(f'Out-of-date export: {target}')
        else:
            target.write_text(expected, encoding='utf-8')
    print('Notebook Python exports match.' if args.check else 'Exported notebook Python cells.')


if __name__ == '__main__':
    main()
