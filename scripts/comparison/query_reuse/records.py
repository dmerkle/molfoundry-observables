"""Verify retained native measurements without depending on repository history.

The snapshot contains the exact measured source bytes. Permitted bookkeeping changes
are relocation of the text formatter import and POSIX manifest-path spelling.
Algorithms and counting code must remain identical.
"""
from hashlib import sha256
from pathlib import Path
import tarfile


def validate_sources(root, summary):
    root = Path(root)
    snapshot = root / "tests/data/measurements/source-snapshot.tar.gz"
    with tarfile.open(snapshot, "r:gz") as archive:
        for name, digest in summary["sources"].items():
            member = archive.extractfile(name)
            if member is None or sha256(member.read()).hexdigest() != digest:
                raise ValueError("Recorded measurement source mismatch: " + name)
            data = (root / name).read_bytes()
            if sha256(data).hexdigest() != digest and name == "src/mfobs/outputs.py":
                data = data.replace(b"from mfobs.formatting import", b"from artifact.show import")
            if name == "src/mfobs/provenance.py":
                data = data.replace(b"p.relative_to(root).as_posix()", b"str(p.relative_to(root))")
            if name == "src/mfobs/core/source_identity.py":
                data = data.replace(b"result[relative.as_posix()]", b"result[str(relative)]")
            if sha256(data).hexdigest() != digest:
                raise ValueError("Measured analysis source changed: " + name)
        script = archive.extractfile("experiments/comparison/query_reuse/measure.py")
        if script is None or sha256(script.read()).hexdigest() != summary["script_sha256"]:
            raise ValueError("Recorded measurement runner mismatch")
    return {"source_commit": summary["source_commit"],
            "verified_sources": len(summary["sources"]),
            "layout_change": "formatter import and POSIX manifest paths; algorithms unchanged"}
