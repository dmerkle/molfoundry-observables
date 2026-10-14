"""Atom-resolved lactone ROP with the ABA query."""
# Example definitions and exact counting on its Ω.
from mfobs.models.molecular_rop.inputs import state_family, configured_variant
from mfobs.models.molecular_rop.inputs import source_rules, query_pattern
from mfobs.models.molecular_rop.counting import MolecularReducer
from mfobs.models.molecular_rop.events import polynomial_jump_model
# Shared framework algorithms.
from mfobs.core.pattern_support import discover_support
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure

MOTIF = "ABA"  # Try "ABBA" with the same chemistry and counting backend.
GATED = False
selection = configured_variant(motifs=(MOTIF,), alphabet=("A", "B"), gated=GATED)
Omega = state_family(selection.alphabet, selection.gated)
rules = source_rules(selection, state_space=Omega)
query = query_pattern(MOTIF)
observables = {query.name: query}

reducer = MolecularReducer(selection, state_space=Omega)
# Algorithm 1: compose the molecular rules and reduce their count images.
linear = discover_linear_observables(rules, observables, counter=reducer, maximum_dimension=24)
# Algorithm 2: discover pattern support and extract polynomial events.
support = discover_support(
    rules, observables, counter=reducer,
)
events = polynomial_jump_model(support, rules, name="molecular-" + MOTIF)
pjc = analyze_polynomial_jump_closure(events)

print("Selected query:", list(observables))
print("Discovered support:", list(support.patterns))
print("Linear closure:", linear["status"])
print("Polynomial jump closure:", pjc["status"])

result = {"linear": linear, "support_discovery": support.record(), "pjc": pjc}
