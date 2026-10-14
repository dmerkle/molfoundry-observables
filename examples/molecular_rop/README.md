# Molecular ROP

The [CLI tutorial](../../tutorials/rop_analysis.py) exposes the motif and gate.
[Inputs](../../src/mfobs/models/molecular_rop/inputs.py),
[counting](../../src/mfobs/models/molecular_rop/counting.py) and
[events](../../src/mfobs/models/molecular_rop/events.py) define the chemical family.

```sh
python -m tutorials rop
python -m mfobs.example_cli rop --out /tmp/rop.json
```
