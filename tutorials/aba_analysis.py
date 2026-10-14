"""ABA counts on a growing path with a switchable append gate."""
# Example definitions and exact counting on its Ω.
from mfobs.models.aba.inputs import state_family, source_rules, query_pattern
from mfobs.models.aba.counting import PathCounter
from mfobs.models.aba.events import polynomial_jump_model
# Shared framework algorithms.
from mfobs.core.pattern_support import discover_support
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure

MOTIF = "ABA"
GATED = False
Omega = state_family(alphabet=("A", "B"), gated=GATED)
rules = source_rules(Omega.definition["alphabet"], Omega.definition["gated"], state_space=Omega)
query = query_pattern(MOTIF)
observables = {query.name: query}

reducer = PathCounter(state_space=Omega)
# Algorithm 1: generator images from rules and exact counts on Ω.
linear = discover_linear_observables(rules, observables, counter=reducer, maximum_dimension=24)

# Algorithm 2: discover pattern support and extract polynomial events.
support = discover_support(
    rules, observables, counter=reducer,
)
events = polynomial_jump_model(support, rules, name="ABA")

pjc = analyze_polynomial_jump_closure(events)

print("Selected query:", list(observables))
print("Discovered support:", list(support.patterns))
print("Linear closure:", linear["status"])
print("Polynomial jump closure:", pjc["status"])

result = {"linear": linear, "support_discovery": support.record(), "pjc": pjc}
