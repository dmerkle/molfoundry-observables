"""Aldol compositions under local chemistry; this example does not synthesize an ODE."""
from collections import Counter
from pathlib import Path

# Framework: native rules, composition and shared chemical restrictions.
import molfoundry as mf
# Example: the local balanced-valence filter and literal GML inputs.
from examples.tcs.inputs import state_family, chemical_rule

Omega = state_family()
folder = Path(__file__).resolve().parents[1] / "examples/tcs"
rules = {name: mf.Rule.fromGMLFile(str(folder / file), add=False) for name, file in
         (("plain", "aldol_plain.gml"), ("noedge", "aldol.gml"))}
observables = {name: mf.Rule.fromGMLFile(str(folder / file), add=False) for name, file in
               (("plain", "query_unconditioned.gml"), ("noedge", "query_conditioned.gml"))}
for rule in (*rules.values(), *observables.values()):
    Omega.validate_rule(rule)

common = mf.rcCommon(connected=False, includeEmpty=True)
evaluator = mf.RCEvaluator([])
results = {}
for name, reaction, query in (
    ("plain", rules["plain"], observables["plain"]),
    ("reaction_noedge", rules["noedge"], observables["plain"]),
    ("reaction_and_query_noedge", rules["noedge"], observables["noedge"]),
):
    after = [r for r in evaluator.eval(reaction * common * query, onlyUnique=False)
             if chemical_rule(r, Omega)]
    before = [r for r in evaluator.eval(query * common * reaction, onlyUnique=False)
              if chemical_rule(r, Omega)]
    net = Counter(r.canonKey for r in after)
    net.subtract(r.canonKey for r in before)
    print(name, "after:", len(after), "before:", len(before),
          "signed terms:", sum(value != 0 for value in net.values()))
    representatives = {rule.canonKey: rule for rule in (*after, *before)}
    results[name] = [{"coefficient": str(coefficient),
                      "gml": representatives[key].getGMLString()}
                     for key, coefficient in net.items() if coefficient]

result = {"compositions": results}
