# Structural NoEdge oscillator

This benchmark is an **artificial but chemically typed** family of covalent
macrocycles.  It is deliberately not an ROP model and does not claim an
experimental mechanism.  Its purpose is to force structural graph reasoning:
Q0/Q1/Q2 are three ordinary C-C closure bonds on one repeated scaffold, and
the guarded channels inspect the absence of one such bond on a neighbouring
site.

The state space is infinite because rings of every size `N >= 7` are admitted,
with `3^N` local topology words before rotational quotienting. Three additional
DPO condensation rules consume a capped repeat, release mapped water, and grow
the ring from `N` to `N+1` without creating atoms inside the graph rule.

`model.json` starts the exact observable worklist from Q0, Q1, Q2.  The reducer
derives connected pair/triplet/longer word patterns from the literal rule
action.  A finite budget must be reported as `UNKNOWN`; the interactive R1-R3
integrations are explicitly labelled Markov approximations.

Rebuild with:

```sh
source /Users/daniel/src/mod-env/bin/activate
python scripts/build_structural_nac_oscillator.py
```
