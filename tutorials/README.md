# Read and run the examples

From the repository root, in the installed analysis environment:

```sh
python -m tutorials
python -m tutorials tutorial
python -m tutorials aba
python -m tutorials rop
python -m tutorials tcs
python -m tutorials chlorination
```

Start with the graph tutorial, then ABA's growing paths. Change `MOTIF` to
`"ABBA"` in [aba_analysis.py](aba_analysis.py). Molecular ROP uses the same
query with atom-resolved rules; TCS demonstrates conditional composition;
chlorination compares a rule with and without its neighbour NAC.

Save a calculation for inspection:

```sh
python -m tutorials eight-sites --out /tmp/eight-sites
```

Start with `summary.txt`: it matches the compact terminal output. `observables.txt`
lists the actual definitions and points to their GML files, including `P33`.
`ode.txt` is present when a finite expectation ODE was established. Query and
rule variants have separate subdirectories. The output directory must be new
or empty; without `--out`, nothing is saved.

For the complete JSON record and available derivation traces, add `--details`:

```sh
python -m tutorials eight-sites --out /tmp/eight-sites-details --details
```

Read the executed code and its actual backend directly:

```sh
python -m tutorials aba --source
python -m tutorials aba --source inputs
python -m tutorials aba --source counting
python -m tutorials aba --source events
```

`mfobs.core` contains the shared algorithms. `mfobs.models.<example>` contains
that family's inputs, exact counting and event arguments. `ABA` and
`MolecularROP` are separate applications; neither chooses a different
representation behind a flag. A new graph family may require its own exact
counting backend. [Adding a model](../docs/new-model.md) explains the boundary.

`python -m tutorials` also lists Rémy, Eight sites, the DNA walker and the
further examples. `comparison` reads retained FRAGGER measurements.
[Notebooks](../notebooks/paper_examples.ipynb) add visualizations;
[the algorithm map](../notebooks/visual/CORE.md) links the shared core.
