"""Symbolic eight-site synthesis from one occupied-slot query."""
from pathlib import Path
from mfobs import ModelInputs, ManifestBackend, analyze
from mfobs.inputs import load

manifest = load(Path(__file__).resolve().parents[1] / "examples/eight-slots/model.json")
defaults = manifest.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
inputs = ModelInputs(Omega, rules, observables)
result = analyze(inputs, backend=ManifestBackend(manifest), max_expansions=100)

print("Selected queries:", list(observables))
print("Discovered basis:", list(result.linear.data["basis"]))
print("Closure:", result.linear.status)
