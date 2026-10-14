# Exact event analysis. Oscillatory macroscopic dynamics in the full notebook are a separate approximation.

# Shared framework input types and dispatcher.
from mfobs import ModelInputs, analyze
# Packaged example and its counting backend.
from mfobs.applications import BrusselatorROP

backend = BrusselatorROP()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
inputs = ModelInputs(Omega, rules, observables)
result = analyze(inputs, backend=backend, linear_limit=24)

print("Selected queries:", list(observables))
print("Linear closure:", result.linear.status)
print("Polynomial jump closure:", result.pjc.status)
