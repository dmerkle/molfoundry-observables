"""Tutorial result saving and protection of existing output directories."""
import csv
import json

import pytest

from tutorials.__main__ import ROOT, main


def test_comparison_saves_its_actual_measurements(tmp_path):
    out = tmp_path / 'comparison'
    assert main(['comparison', '--out', str(out)]) == 0
    with (ROOT / 'notebooks/visual/evidence/ab-scaling/summary.csv').open(encoding="utf-8") as stream:
        expected = list(csv.DictReader(stream))
    with (out / 'measurements.csv').open(encoding="utf-8") as stream:
        assert list(csv.DictReader(stream)) == expected
    assert 'FRAGGER seconds' in (out / 'summary.txt').read_text(encoding="utf-8")
    assert not (out / 'result.json').exists()
    assert not (out / 'ode.txt').exists()


def test_brusselator_counting_source_is_the_executed_backend(capsys):
    import inspect
    from pathlib import Path
    from mfobs.models.brusselator_growth_projection import BrusselatorGrowthGraphModel

    assert main(['brusselator', '--source', 'counting']) == 0
    actual = Path(inspect.getsourcefile(BrusselatorGrowthGraphModel))
    assert capsys.readouterr().out == actual.read_text(encoding="utf-8")
    assert actual.name in (ROOT / 'notebooks/visual/CORE.md').read_text(encoding="utf-8")


def test_full_records_require_details(tmp_path):
    out = tmp_path / 'comparison'
    assert main(['comparison', '--out', str(out), '--details']) == 0
    with (ROOT / 'notebooks/visual/evidence/ab-scaling/summary.csv').open(encoding="utf-8") as stream:
        expected = list(csv.DictReader(stream))
    assert json.loads((out / 'result.json').read_text(encoding="utf-8")) == {'measurements': expected}


def test_aba_compact_export_retains_patterns_but_claims_no_closed_ode(tmp_path, aba_records):
    from mfobs.outputs import export_summary
    record = aba_records['aba']
    export_summary(record, tmp_path, name='ABA')
    summary = (tmp_path / 'summary.txt').read_text(encoding="utf-8")
    assert 'No finite closed linear ODE established' in summary
    assert 'Event classes: 8' in summary
    linear, events = summary.split('Finite events and EMGF (Algorithm 2)')
    assert 'd⟨' not in linear
    assert 'Exact first-moment identities from these events' in events
    assert len(summary.splitlines()) < 70
    assert not (tmp_path / 'result.json').exists()
    assert not (tmp_path / 'trace.json').exists()
    assert not (tmp_path / 'ode.txt').exists()
    index = (tmp_path / 'observables.txt').read_text(encoding="utf-8")
    assert 'C_ABA = (1) count(observables/C_ABA-0.gml)' in index
    assert 'M_A*S_AB = M_A*S_AB' not in index
    assert (tmp_path / 'observables/C_ABA-0.gml').read_text(encoding="utf-8") == record['support_discovery']['patterns']['C_ABA']['gml']


def test_eight_site_computed_definitions_and_ode_are_immediately_accessible(tmp_path, capsys):
    out = tmp_path / 'eight-sites'
    assert main(['eight-sites', '--out', str(out)]) == 0
    terminal = capsys.readouterr().out
    assert 'Status: closed modulo eight slot axioms' in terminal
    assert '48 retained observables' in terminal
    assert 'Statistics: compositions:' in terminal
    assert len(terminal.splitlines()) < 30
    assert not (out / 'result.json').exists()
    assert not (out / 'trace.json').exists()
    index = (out / 'observables.txt').read_text(encoding="utf-8")
    assert 'P33 = ' in index and 'count(observables/P33-0.gml)' in index
    assert 'constrainNoEdge' in (out / 'observables/P33-0.gml').read_text(encoding="utf-8")
    assert 'dmu_P33/dt =' in (out / 'ode.txt').read_text(encoding="utf-8")


def test_all_query_variation_results_have_their_own_definitions(tmp_path, aba_records):
    from mfobs.outputs import export_summary
    aba = aba_records['aba']
    export_summary({'queries': {'ABA': aba, 'AB': aba}}, tmp_path)
    assert '  ABA |' in (tmp_path / 'summary.txt').read_text(encoding="utf-8")
    assert '  AB |' in (tmp_path / 'summary.txt').read_text(encoding="utf-8")
    assert len((tmp_path / 'summary.txt').read_text(encoding="utf-8").splitlines()) == 4
    for name in ('ABA', 'AB'):
        assert (tmp_path / name / 'observables/C_ABA-0.gml').exists()


def test_native_dna_walker_retains_its_derived_ode_and_sparse_definitions(tmp_path, capsys):
    out = tmp_path / 'walker'
    assert main(['dna-walker', '--out', str(out)]) == 0
    text = capsys.readouterr().out
    assert 'Selected queries: E' in text
    assert 'Exact finite ODE' in text
    assert 'd⟨E⟩/dt' in (out / 'ode.txt').read_text(encoding="utf-8")
    assert 'E = (1) count(observables/E-0.gml)' in (out / 'observables.txt').read_text(encoding="utf-8")
    assert (out / 'observables/E-0.gml').exists()


def test_output_directory_is_checked_before_execution(tmp_path, monkeypatch):
    old = tmp_path / 'result.json'
    old.write_text('keep this calculation', encoding="utf-8")
    def unexpected_run(*args, **kwargs):
        pytest.fail('an invalid output directory must be rejected before calculation')
    monkeypatch.setattr('tutorials.__main__.runpy.run_module', unexpected_run)
    with pytest.raises(SystemExit) as error:
        main(['eight-sites', '--out', str(tmp_path)])
    assert error.value.code == 2
    assert old.read_text(encoding="utf-8") == 'keep this calculation'


@pytest.mark.parametrize('arguments', [
    ['--out', 'unused'],
    ['eight-sites', '--source', '--out', 'unused'],
    ['aba', '--details'],
    ['--details'],
])
def test_output_requires_execution(arguments):
    with pytest.raises(SystemExit) as error:
        main(arguments)
    assert error.value.code == 2


def test_support_export_does_not_leave_an_earlier_ode(tmp_path):
    from mfobs.outputs import export
    for name in ('ode.txt', 'ode.json', 'matrix.csv', 'trace.json'):
        (tmp_path / name).write_text('old linear system', encoding="utf-8")
    patterns = tmp_path / 'observables'
    patterns.mkdir()
    (patterns / 'old.gml').write_text('old pattern', encoding="utf-8")
    record = {'pjc': {'status': 'UNKNOWN'}, 'support_discovery': {'patterns': {}}}
    export(record, tmp_path)
    assert json.loads((tmp_path / 'result.json').read_text(encoding="utf-8")) == record
    assert {path.name for path in tmp_path.iterdir()} == {'result.json'}
