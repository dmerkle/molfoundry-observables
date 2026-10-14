#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON=${PYTHON:-python}
"$PYTHON" run.py doctor
"$PYTHON" run.py run examples/gated-path/model.json --out out/gated-path --force --progress
"$PYTHON" run.py verify examples/gated-path/model.json --result out/gated-path/result.json --full --multiplicities --conditions
"$PYTHON" run.py run examples/ring8/model.json --out out/ring8 --force --progress
"$PYTHON" run.py verify examples/ring8/model.json --result out/ring8/result.json --multiplicities
"$PYTHON" -m pytest -q
