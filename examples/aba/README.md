# Growing-path query

The [CLI tutorial](../../tutorials/aba_analysis.py) defines Ω, source rules and
the selected motif. [Inputs](../../src/mfobs/models/aba/inputs.py),
[counting](../../src/mfobs/models/aba/counting.py) and
[events](../../src/mfobs/models/aba/events.py) contain the actual family code.

`query-gml.json` selects the literal pattern in `observables/`; `query.json`
contains the default query and search budgets for JSON exports.

```sh
python -m tutorials aba
python -m mfobs.example_cli aba --out /tmp/aba.json
```
