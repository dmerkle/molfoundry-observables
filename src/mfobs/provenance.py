"""Fingerprints for checking whether computed records can be reused."""
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, distribution, version
from pathlib import Path
import json
import platform


def source_hashes():
    root = Path(__file__).resolve().parents[2]
    presentation = {'chemical_depiction.py', 'rule_drawing.py', 'visual_links.py', 'cli.py', '__main__.py'}
    paths = [p for p in (root / 'src/mfobs').rglob('*.py') if p.name not in presentation]
    paths.extend((root / 'src/mfobs/data').glob('*.json'))
    paths.append(root / 'pyproject.toml')
    return {p.relative_to(root).as_posix(): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def runtime_fingerprint():
    def installed(name):
        try:
            return version(name)
        except PackageNotFoundError:
            return None
    dist = distribution('molfoundry')
    paths = sorted(p for p in (dist.files or ())
                   if str(p).startswith(('molfoundry/', '_molfoundry_pm/'))
                   and p.suffix in {'.py', '.js', '.wasm', '.so', '.dylib', '.dll'})
    files = {str(p): sha256(Path(dist.locate_file(p)).read_bytes()).hexdigest() for p in paths}
    return {
        'python': platform.python_version(),
        'versions': {name: installed(name) for name in ('molfoundry', 'sympy', 'networkx', 'numpy', 'scipy')},
        'molfoundry_sha256': sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
    }


def check_runtime(record):
    if record.get('runtime') != runtime_fingerprint():
        raise ValueError('Analysis runtime changed or is unrecorded; recompute the results')
