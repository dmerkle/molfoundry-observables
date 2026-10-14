# Retained FRAGGER comparison

`walker-export.json` contains the original DNA-walker rules, query and FRAGGER
equations. `walker-incidence.json` contains the incidence-encoded comparison.
The two Scala harnesses are the executed source. `../result.json` and
`../provenance.json` bind the source and result hashes.

Retained from development CI run `35866705340`,
source commit `ab432a9d6c48d05d6bc758621bb370aaedf660d0`. The current notebook
recomputes MolFoundry's equations; it reads these files only for comparison.
