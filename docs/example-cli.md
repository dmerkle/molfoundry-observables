# Inspect inputs and run analyses

The [package README](../README.md) explains installation, Docker and notebooks. The [algorithm map](../notebooks/visual/CORE.md) points to the shared residual loop, support discovery and backend-specific proof obligations.

The CLI dispatcher is [mfobs.example_cli](../src/mfobs/example_cli.py).
Literal inputs live in `examples/`;
Python model definitions and counting backends live in `src/mfobs/models/`.
From the package root:

```sh
PYTHONPATH=src python -m mfobs.example_cli aba --out /tmp/aba.json
PYTHONPATH=src python -m mfobs.example_cli remy --out /tmp/remy.json
PYTHONPATH=src python -m mfobs.example_cli chlorination-guarded --out /tmp/chlorination.json
```

To vary a supported query, use `--query-file` or repeat `--observable NAME=FILE[@COEFFICIENT]`. Repeated names form one finite rule expression; Rémy's P3 query includes three stem variants. The result retains selected files and hashes under `query_input`. Unsupported query shapes must be rejected rather than reinterpreted.

```sh
PYTHONPATH=src python -m mfobs.example_cli remy \
  --query-file examples/remy/query.json --out /tmp/remy-query.json
```

For the paper's symbolic eight-site analysis use [examples/eight-slots](../examples/eight-slots/README.md). The `ring8` command instead evaluates a fixed finite quotient. Aldol is a conditional-composition experiment, not an ODE-closure example. The [DNA-walker notebook](../notebooks/examples/06_dna_walker.ipynb) recomputes MolFoundry's discovery through the shared linear-closure algorithm and compares it with retained FRAGGER results.

Model construction and reduction assumptions are inputs. Names assigned to discovered coordinates are presentation; inspect their literal expressions and discovery origins. Event extraction may require an additional uniform-match argument or finite-language construction even when support discovery succeeds. Budget exhaustion is reported as UNKNOWN, not non-closure.
