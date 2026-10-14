"""Keep the example excerpts tied to runnable, inspectable examples."""
import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_standalone_catalogue_matches_visual_notebooks():
    examples = json.loads((ROOT / 'tutorials/minimal_examples.json').read_text(encoding='utf-8'))
    active = ROOT / 'notebooks/examples'
    assert {row['notebook'] for row in examples if row['notebook'] is not None} == {
        p.stem for p in active.glob('[0-9][0-9]_*.ipynb')
    }
    for row in examples:
        script = ROOT / 'tutorials' / (row['script'] + '.py')
        code = script.read_text(encoding='utf-8').strip()
        ast.parse(code)
        if row['notebook'] is None:
            continue
        notebook = json.loads((active / (row['notebook'] + '.ipynb')).read_text(encoding='utf-8'))
        assert not any(
            '<!-- minimal-inspection -->' in ''.join(cell['source'])
            for cell in notebook['cells']
        )  # Standalone scripts are not duplicated in the visual walkthrough.


def test_minimal_examples_have_no_plotting_or_notebook_execution_dependency():
    examples = json.loads((ROOT / 'tutorials/minimal_examples.json').read_text(encoding='utf-8'))
    forbidden = {'matplotlib', 'ipywidgets', 'paper_view', 'visual_view', 'visual_suite', 'nbclient'}
    for row in examples:
        tree = ast.parse((ROOT / 'tutorials' / (row['script'] + '.py')).read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or '']
            assert not any(set(name.split('.')) & forbidden for name in names), row['script']
