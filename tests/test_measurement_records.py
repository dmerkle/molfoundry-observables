"""Integrity of the compact native measurements shipped with the source."""
import gzip
import hashlib
import json
from pathlib import Path
import statistics
import shutil
import pytest
from scripts.comparison.query_reuse.records import validate_sources

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "tests/data/measurements"


def test_query_records_bind_to_current_source_and_complete_raw_results():
    summary = json.loads((RECORDS / "query-discovery.json").read_text(encoding="utf-8"))
    assert summary["status"] == "PASS" and len(summary["observations"]) == 15
    assert validate_sources(ROOT, summary)["verified_sources"] == len(summary["sources"])
    for row in summary["observations"]:
        data = gzip.decompress((RECORDS / "query-runs" /
            f"{row['motif']}-{row['repetition']}.json.gz").read_bytes())
        assert hashlib.sha256(data).hexdigest() == row["raw_sha256"]
        assert row["oracle"]["status"] == "PASS"
    for motif, record in summary["cases"].items():
        rows = [row for row in summary["observations"] if row["motif"] == motif]
        assert len(rows) == 3
        assert record["discovery_and_events_seconds"]["median"] == statistics.median(
            row["discovery_and_events_seconds"] for row in rows)


def test_fragger_measurements_are_retained_unchanged():
    comparison = json.loads((RECORDS / "composition.json").read_text(encoding="utf-8"))
    original = (RECORDS / "fragger-original.json").read_bytes()
    assert hashlib.sha256(original).hexdigest() == comparison["retained_fragger"]["sha256"]
    rows = lambda data: [row for row in data["observations"] if row["tool"] == "fragger-core"]
    assert rows(comparison) == rows(json.loads(original))


@pytest.mark.parametrize("source", ["core/worklist.py", "core/source_identity.py", "provenance.py"])
def test_recorded_measurements_reject_changed_algorithm_source(tmp_path, source):
    summary = json.loads((RECORDS / "query-discovery.json").read_text(encoding="utf-8"))
    for name in summary["sources"]:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
    target = tmp_path / "tests/data/measurements/source-snapshot.tar.gz"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(RECORDS / "source-snapshot.tar.gz", target)
    algorithm = tmp_path / "src/mfobs" / source
    algorithm.write_bytes(algorithm.read_bytes() + b"\n# source drift\n")
    with pytest.raises(ValueError, match="Measured analysis source changed"):
        validate_sources(tmp_path, summary)
