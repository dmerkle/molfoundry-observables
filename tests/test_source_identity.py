import json
from pathlib import Path

import pytest

from mfobs.core.source_identity import (
    make_snapshot_manifest,
    require_source_identity,
    source_identity,
)


def frozen_tree(tmp_path: Path) -> Path:
    root = tmp_path / "source"
    root.mkdir()
    (root / "a.py").write_text("answer = 42\n", encoding="utf-8")
    sub = root / "data"
    sub.mkdir()
    (sub / "model.json").write_text('{"model": true}\n', encoding="utf-8")
    manifest = make_snapshot_manifest(root, source_commit="abc123")
    assert set(manifest["files"]) == {"a.py", "data/model.json"}
    (root / ".mfobs-snapshot.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root


def test_standalone_snapshot_verifies_without_git(tmp_path):
    identity = source_identity(frozen_tree(tmp_path))
    assert identity["mode"] == "verified-frozen-snapshot"
    assert identity["commit"] == "abc123"
    assert identity["verified"] is True
    assert identity["file_count"] == 2


@pytest.mark.parametrize("mutation", ["changed", "missing", "extra"])
def test_snapshot_fails_closed_for_tampered_missing_or_extra_file(tmp_path, mutation):
    root = frozen_tree(tmp_path)
    if mutation == "changed":
        (root / "a.py").write_text("answer = 0\n", encoding="utf-8")
    elif mutation == "missing":
        (root / "a.py").unlink()
    else:
        (root / "extra.txt").write_text("not inventoried\n", encoding="utf-8")
    with pytest.raises(ValueError, match="verification failed"):
        source_identity(root)


def test_missing_manifest_is_not_treated_as_a_clean_snapshot(tmp_path):
    root = tmp_path / "plain"
    root.mkdir()
    (root / "file.py").write_text("pass\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="neither a Git checkout"):
        source_identity(root)


def test_expected_commit_is_enforced(tmp_path):
    root = frozen_tree(tmp_path)
    with pytest.raises(ValueError, match="source commit mismatch"):
        require_source_identity(root, expected_commit="other")
