"""Rémy's six insertion rules and one selected left-spine count."""

# Example definitions, count identities and event argument.
import mfobs.models.remy.counting as remy_counting
import mfobs.models.remy.events as remy_events
import molfoundry as mf

from mfobs.models.remy.inputs import state_family
from mfobs.models.remy import inputs as remy
# Shared framework algorithms.
from mfobs.core.jump_closure_worklist import RuleExpression, discover_observables
from mfobs.core.typed_incidence import pair_commutator
from mfobs.core.polynomial_jump_closure import analyze_polynomial_jump_closure

Omega = state_family()
rules = remy.insertion_rules()
query = tuple(mf.Rule.fromGMLString(remy.observable_gml("P3", stem), add=False)
              for stem in remy.EDGE_TYPES)
observables = {"P3": RuleExpression.from_rules(query)}

linear = discover_observables(
    seeds=observables, generator=RuleExpression.from_rules(rules.values()),
    pair_commutator=pair_commutator(Omega.definition["incidence"], include_empty=True),
    classify=remy_counting.classify_observable,
    jump_closure=remy_counting.valid_tree_dpo_jump_closure,
    maximum_dimension=16,
)
# Complete events additionally use the all-tree argument for these six rules.
events = remy_events.polynomial_jump_model(3)
pjc = analyze_polynomial_jump_closure(events)

print("Supplied rules:", list(rules))
print("Selected query:", list(observables))
print("Discovered basis:", list(linear.basis))
print("Linear closure:", linear.status)
print("Polynomial jump closure:", pjc["status"])

result = {"linear": linear.record(), "pjc": pjc,
          "observable_expressions": {name: [
              {"coefficient": str(coefficient), "gml": expression.representatives[key].getGMLString()}
              for key, coefficient in expression.coefficients.items()
          ] for name, expression in linear.basis.items()}}
