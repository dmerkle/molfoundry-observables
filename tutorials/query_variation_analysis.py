"""Five motif queries with the same ungated path rules and counting backend."""
# Example definitions and exact counting on its Ω.
from mfobs.models.aba.inputs import state_family, source_rules, query_pattern
from mfobs.models.aba.counting import PathCounter
from mfobs.models.aba.events import polynomial_jump_model
# Shared framework algorithms.
from mfobs.core.pattern_support import discover_support
from mfobs.core.counting_linear_closure import discover_linear_observables
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure

Omega = state_family(alphabet=("A", "B"), gated=False)
rules = source_rules(Omega.definition["alphabet"], state_space=Omega)
reducer = PathCounter(state_space=Omega)

results = {}
for motif in ("AB", "ABA", "ABBA", "AAAA", "ABABA"):
    query = query_pattern(motif)
    observables = {query.name: query}
    # Algorithm 1, independently of polynomial support and events.
    linear = discover_linear_observables(rules, observables, counter=reducer, maximum_dimension=24)
    # Algorithm 2.
    support = discover_support(
        rules, observables, counter=reducer,
    )
    events = polynomial_jump_model(support, rules, name=motif)
    pjc = analyze_polynomial_jump_closure(events)
    print(motif, "support:", list(support.patterns))
    print("Linear:", linear["status"], "pJC:", pjc["status"])
    results[motif] = {"linear": linear, "support_discovery": support.record(), "pjc": pjc}

result = {"queries": results}
