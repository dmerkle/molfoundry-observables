"""Semantic package data participates in saved-result validation."""
from mfobs import provenance


def test_source_hashes_include_changes_to_chemical_theory(tmp_path, monkeypatch):
    source = tmp_path / "src/mfobs"
    data = source / "data"
    data.mkdir(parents=True)
    (tmp_path / "pyproject.toml").write_text("[project]\n", encoding="utf-8")
    (source / "provenance.py").write_text("# fixture\n", encoding="utf-8")
    theory = data / "chemistry.json"
    theory.write_text('{"maxValence": {"C": 4}}', encoding="utf-8")
    monkeypatch.setattr(provenance, "__file__", str(source / "provenance.py"))
    before = provenance.source_hashes()
    theory.write_text('{"maxValence": {"C": 3}}', encoding="utf-8")
    after = provenance.source_hashes()
    key = "src/mfobs/data/chemistry.json"
    assert key in before
    assert before[key] != after[key]
    assert before["src/mfobs/provenance.py"] == after["src/mfobs/provenance.py"]
