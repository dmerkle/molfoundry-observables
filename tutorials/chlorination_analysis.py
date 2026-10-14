"""Single-site chlorination, with and without the neighbour NAC."""
# Shared framework input types and dispatcher.
from mfobs import ModelInputs, Observable, RatedRule, analyze
# Example definitions and its specialised counting/event backend.
from mfobs.models.chlorination.analysis import Chlorination
from mfobs.models.chlorination.inputs import chlorination_rule, query_patterns, state_family, grow_rule, RATE_FACTORS

Omega = state_family(chlorine_total=4)
observables = {"Z": Observable("Z", query_patterns()["Z"].rule)}

results = {}
for variant in ("guarded", "ablated"):
    reaction = chlorination_rule(guarded=variant == "guarded")
    # 3! terminal-H permutations × 4! methane-H permutations: rate g per molecule pair.
    rules = {"chlorinate": RatedRule(reaction, "k", identifier="chlorinate"),
             "grow": RatedRule(grow_rule(), "g", per_match=RATE_FACTORS["grow"])}
    inputs = ModelInputs(Omega, rules, observables)
    result = analyze(inputs, backend=Chlorination(variant))
    print(variant, result.pjc.status)
    print("Discovered counts:", result.pjc.data["coordinates"])
    results[variant] = result.record()

result = {"variants": results}
