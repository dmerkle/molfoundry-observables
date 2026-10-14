# Bounded exact search on an unbounded ring language. A budget stop is not a non-closure proof; R1 dynamics are a separate approximation.

# Shared framework input types and dispatcher.
from mfobs import ModelInputs, analyze
# Packaged example and its counting backend.
from mfobs.applications import GrowingRing

backend = GrowingRing()
defaults = backend.specification()
Omega = defaults.omega
rules = defaults.rules
observables = defaults.observables
inputs = ModelInputs(Omega, rules, observables)
# Same bounded prefix as the notebook; increase for a longer exact search.
result = analyze(inputs, backend=backend, worklist_limit=12)

print("Selected queries:", list(observables))
print("Linear closure:", result.linear.status)
print("Polynomial jump closure:", result.pjc.status)
print("Scope:", result.pjc.data["reason"])
