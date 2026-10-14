# Rémy trees: six insertion rules, one expression-valued seed

[`query.json`](query.json) is the same common query format used by the other
examples. It supplies **one** observable, `P3`, as the sum of the three literal
identity patterns [`P3-I.gml`](observables/P3-I.gml),
[`P3-L.gml`](observables/P3-L.gml), and
[`P3-R.gml`](observables/P3-R.gml). The three files distinguish whether the
top occurrence is attached through the exterior stem, a left edge, or a right
edge; together they count every occurrence once.

[`run.py`](run.py) calls the two-branch analysis in
[`remy/analysis.py`](../../src/mfobs/models/remy/analysis.py).
[`remy/inputs.py`](../../src/mfobs/models/remy/inputs.py) defines Ω and generates
the six GML insertion rules from `side ∈ {L,R}` and incoming edge type `I,L,R`.
[`remy/counting.py`](../../src/mfobs/models/remy/counting.py) supplies the exact
tree identities. Linear discovery starts from **P3 only**;
the EMGF question explicitly selects all four counts. There is no saved
jump table passed to the generic pJC compiler.

The generic worklist compares canonical sparse rule expressions. It neither
needs nor uses the names `P2`, `P1`, or `E` to decide equality, independence,
or closure. Only after a new residual has been accepted does the small
tree-language recognizer attach the conventional `P_k` interpretation; an
uninterpreted generic name would yield the same span. Read
`observable_designation` for this boundary.

Read `contract` for exact rule and observable declarations, `composition`
for the scope of rule-composition evidence, `linear` for the discovered
basis, and `pjc` for event classes and operator. The all-tree activity
identity and the P2/P3 context partition are model-specific combinatorial
arguments, not full generic graph-composition output. Their exact statement,
bijection and proof boundary are retained in `model_specific_event_reduction`.
The artifact profile
also runs independent bounded recursive-tree and direct MolFoundry checks;
these validate but do not replace the all-tree argument.

```sh
PYTHONPATH=src python -m mfobs.example_cli remy --out /tmp/remy.json
PYTHONPATH=src python -m mfobs.example_cli remy \
  --query-file examples/remy/query.json --out /tmp/remy-from-gml.json
```
