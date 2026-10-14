#!/usr/bin/env python3
"""Define Omega, one GML rule and one query; run Algorithms 1 and 2 separately.

Run from:

    # From the repository root:
    PYTHONPATH=src python tutorials/graph_closure_tutorial.py
"""

from molfoundry import ruleGMLString, graphGMLString
from mfobs import ModelInputs, RatedRule, Observable, SimpleGraphs, SimpleGraphBackend


# 1) State language Omega:
#    all finite simple graphs whose vertices are L/R and whose edges are "-".
Omega = SimpleGraphs(
    vertex_labels={"L", "R"},
    edge_labels={"-"},
)


# 2) One DPO graph rule:
#
#       L -- R     --->     L -- L -- R
#
subdivide = ruleGMLString("""rule [
  ruleID "subdivide"
  context [
    node [ id 0 label "L" ]
    node [ id 1 label "R" ]
  ]
  left [ edge [ source 0 target 1 label "-" ] ]
  right [
    node [ id 2 label "L" ]
    edge [ source 0 target 2 label "-" ]
    edge [ source 2 target 1 label "-" ]
  ]
]""")
rules = {"subdivide": RatedRule(subdivide, rate="k")}


# 3) The supplied observable:
#
#       N_L = number of matches of the one-vertex graph [L]
#
N_L = Observable(
    "N_L",
    graphGMLString('graph [ node [ id 0 label "L" ] ]'),
)


observables = {N_L.name: N_L}
inputs = ModelInputs(Omega, rules, observables)
backend = SimpleGraphBackend(Omega)


# 4) Algorithm 1: discover a closed linear observable system.
linear = backend.linear_closure(inputs)
print("Algorithm 1:", linear.status)
print("Basis:", list(linear.basis))
for step in linear.steps:
    print(f"G[{step['channel']}] {step['source']} =", step["coordinates"])


# 5) Algorithm 2: discover pattern support and polynomial jump events.
# Uniformity: every L--R subdivision creates one L and preserves the L--R count.
jump = backend.jump_closure(inputs, uniform_matches=True)
print("\nAlgorithm 2:", jump.pjc["status"])
print("Support:", list(jump.support.patterns))
if jump.compilation is not None:
    for event in jump.compilation.model.events:
        activity = jump.compilation.model.ring.text(event.propensity)
        print(f"Activity: {event.channel}*({activity}); jump: {event.delta}")

# Both discover E_LR = count([L--R]):
# G N_L = k E_LR, G E_LR = 0; event activity k E_LR and jump (+1, 0).

result = {"linear": linear.record(), "support_discovery": jump.support.record(), "pjc": jump.pjc}
