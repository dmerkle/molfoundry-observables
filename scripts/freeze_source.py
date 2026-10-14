"""Record the exact source copied into a Docker image (not an analysis cache)."""
from pathlib import Path
import json
import sys
from mfobs.core.source_identity import make_snapshot_manifest, verify_snapshot

root = Path(__file__).resolve().parents[1]
manifest = make_snapshot_manifest(
    root, source_commit=sys.argv[1],
    ignored_directories=('molfoundry_observables.egg-info', '.pytest_cache',
                         '.ipynb_checkpoints', 'runs', 'inspection'))
(root / '.mfobs-snapshot.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding="utf-8")
verify_snapshot(root)
