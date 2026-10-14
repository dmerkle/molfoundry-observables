"""Presentation/optional diagnostics must not block or misrepresent synthesis."""
from pathlib import Path
import json
import shutil
from mfobs.inputs import load
from mfobs.core.worklist import run

ROOT = Path(__file__).resolve().parents[1]

def test_krylov_determinant_can_be_disabled_without_changing_closure(tmp_path):
    target=tmp_path/'model'
    shutil.copytree(ROOT/'examples/gated-path',target)
    path=target/'reduction.json'
    data=json.loads(path.read_text(encoding='utf-8'));data['krylov_max_dimension']=0
    path.write_text(json.dumps(data), encoding='utf-8')
    result=run(load(target/'model.json'))
    assert result['closed'] and len(result['basis'])==3
    assert 'krylov_determinant' not in result
    assert result['krylov_diagnostic']['status']=='skipped'


def test_verification_identifies_declared_grammar_hosts():
    from mfobs.checks.invariant_discovery import verify
    inputs=load(ROOT/'examples/gated-path/model.json')
    result=run(inputs)
    checked=verify(inputs,result)
    assert checked['passed'] and checked['host_constructor']=='declared-grammar'
    assert checked['independent_host_constructor'] is False
