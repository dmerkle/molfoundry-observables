"""Honest source identity for development checkouts and frozen snapshots."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Iterable


PROJECT_SNAPSHOT = ".mfobs-snapshot.json"
FRAGGER_SNAPSHOT = ".fragger-snapshot.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _files(root: Path, *, snapshot_name: str,
           ignored_directories: Iterable[str] = ()) -> dict[str, str]:
    ignored = set(ignored_directories)
    result = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if (not path.is_file() or relative.name == snapshot_name
                or "__pycache__" in relative.parts or path.suffix == ".pyc"
                or any(part in ignored for part in relative.parts)):
            continue
        result[relative.as_posix()] = _sha256(path)
    return result


def make_snapshot_manifest(root: Path, *, source_commit: str,
                           snapshot_name: str = PROJECT_SNAPSHOT,
                           ignored_directories: Iterable[str] = ()) -> dict:
    """Create a complete manifest for an already staged source directory."""

    files = _files(root, snapshot_name=snapshot_name,
                   ignored_directories=ignored_directories)
    return {
        "schema_version": 1,
        "kind": "verified-frozen-source-snapshot",
        "source_commit": source_commit,
        "snapshot_file": snapshot_name,
        "ignored_generated_directories": list(ignored_directories),
        "files": files,
        "inventory_sha256": _digest(files),
    }


def verify_snapshot(root: Path, *, snapshot_name: str = PROJECT_SNAPSHOT) -> dict:
    path = root / snapshot_name
    if not path.is_file():
        raise FileNotFoundError(f"frozen source manifest is missing: {path}")
    record = json.loads(path.read_text(encoding="utf-8"))
    if (record.get("schema_version") != 1
            or record.get("kind") != "verified-frozen-source-snapshot"
            or record.get("snapshot_file") != snapshot_name
            or not record.get("source_commit")):
        raise ValueError("invalid frozen source manifest schema")
    expected = record.get("files")
    if not isinstance(expected, dict) or not expected:
        raise ValueError("frozen source manifest has no file inventory")
    if record.get("inventory_sha256") != _digest(expected):
        raise ValueError("frozen source inventory digest is invalid")
    actual = _files(
        root, snapshot_name=snapshot_name,
        ignored_directories=record.get("ignored_generated_directories", ()),
    )
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    changed = sorted(name for name in set(expected) & set(actual)
                     if expected[name] != actual[name])
    if missing or extra or changed:
        raise ValueError(
            "frozen source verification failed: "
            f"missing={missing}, extra={extra}, changed={changed}"
        )
    return {
        "mode": "verified-frozen-snapshot",
        "commit": record["source_commit"],
        "verified": True,
        "snapshot_manifest": str(path),
        "inventory_sha256": record["inventory_sha256"],
        "file_count": len(actual),
        "ignored_generated_directories": record.get("ignored_generated_directories", []),
    }


def git_identity(root: Path, *, paths: Iterable[str] = ()) -> dict:
    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments], cwd=root, text=True, capture_output=True, check=True,
        ).stdout.strip()

    scope = tuple(paths)
    status = git("status", "--porcelain", *(('--', *scope) if scope else ()))
    return {
        "mode": "git-checkout",
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_porcelain": status.splitlines(),
        "verified": not bool(status),
        "scope": list(scope) if scope else ["complete checkout"],
    }


def source_identity(root: Path, *, snapshot_name: str = PROJECT_SNAPSHOT,
                    git_paths: Iterable[str] = ()) -> dict:
    """Use an actual Git identity, or verify a complete frozen inventory."""

    root = root.resolve()
    snapshot = root / snapshot_name
    if snapshot.is_file():
        return verify_snapshot(root, snapshot_name=snapshot_name)
    try:
        return git_identity(root, paths=git_paths)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            f"{root} is neither a Git checkout nor a verified frozen snapshot"
        ) from exc


def require_source_identity(root: Path, *, expected_commit: str | None = None,
                            snapshot_name: str = PROJECT_SNAPSHOT,
                            git_paths: Iterable[str] = ()) -> dict:
    identity = source_identity(root, snapshot_name=snapshot_name, git_paths=git_paths)
    if expected_commit is not None and identity["commit"] != expected_commit:
        raise ValueError(
            f"source commit mismatch: {identity['commit']} != {expected_commit}"
        )
    if not identity["verified"]:
        raise ValueError("source identity is dirty or unverified")
    return identity
