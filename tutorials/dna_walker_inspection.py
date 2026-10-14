"""DNA-walker closure from four rules and one query."""
# Shared framework input types and dispatcher.
from mfobs import ModelInputs, analyze
# Packaged example and its counting backend.
from mfobs.applications import DNAWalker

backend = DNAWalker()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
inputs = ModelInputs(Omega, rules, observables)
result = analyze(inputs, backend=backend)

print("Selected query:", list(observables))
print("Discovered basis:", list(result.linear.data["basis"]))
print("Linear closure:", result.linear.status)
print("Normalized equations:", result.evidence["comparison"]["physical_equations"])
