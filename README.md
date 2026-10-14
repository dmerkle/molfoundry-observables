# Observable analysis with MolFoundry

MolFoundry is a graph transformation framework. Its Python analysis package
`mfobs` derives observable dynamics from graph rules, rates and query patterns,
relative to a declared graph family and exact counting backend.

This repository accompanies the paper:

Marcel Friedrichs, Addie Jordon, Juri Kolčák, and Daniel Merkle.
*Automating Observable Closure Analysis for Graph Transformation Models with MolFoundry*.
arXiv link forthcoming.

**Try it online:** [Interactive website and executable Jupyter notebooks](https://stochmech.algochem.techfak.de/).
Explore worked examples on the website or run their analyses in the hosted notebooks.

## Start here

1. Read the [graph tutorial](tutorials/graph_closure_tutorial.py): graph family,
   rule and query are inline; Algorithms 1 and 2 are called independently.
2. Run [ABA](tutorials/aba_analysis.py), then [molecular ROP](tutorials/rop_analysis.py).
3. Change ABA's `MOTIF` from `"ABA"` to `"ABBA"` and rerun.

```sh
python -m tutorials tutorial
python -m tutorials aba
python -m tutorials aba --source counting
python -m tutorials eight-sites --out /tmp/eight-sites
```

`--out` writes `summary.txt`, `observables.txt`, literal GML and, when established,
`ode.txt`. Use a new or empty directory. `--details` adds the full JSON record.
[CLI guide](tutorials/README.md) · [Algorithms and source map](notebooks/visual/CORE.md)

`src/mfobs/core` contains shared algorithms; `src/mfobs/models/<example>` contains
family-specific inputs, counting and event arguments. [Adding a model](docs/new-model.md)
explains that boundary. Budget exhaustion reports UNKNOWN, not non-closure.

## Install and open

Use a cloned checkout: the examples and notebooks are source-tree resources.
Python 3.13 on Linux x86-64 is the pinned Docker environment; compatible local
MolFoundry installations can also run the source.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --extra-index-url https://py.molfoundry.de/dl/ -r requirements.txt -e .
python -m ipykernel install --user --name molfoundry --display-name MolFoundry
jupyter lab notebooks/paper_examples.ipynb
```

In the index, open **02_aba.ipynb**. Choose the **MolFoundry** kernel, then
**Restart Kernel and Run All**.
All eleven numbered notebooks are independent. GitHub shows code and saved
results; it does not execute them. RDKit supplies bond/H highlights and charged
aromatic views; ordinary graph and molecule drawings use MolFoundry.

## Docker

Build from this directory. The build downloads pinned dependencies; the image
then runs the included examples without network access.

```sh
docker build --platform linux/amd64 --build-arg SOURCE_COMMIT="$(git rev-parse HEAD)" -t molfoundry-observables .
docker run --rm --network none molfoundry-observables aba
mkdir -p output
docker run --rm --network none -v "$PWD/output:/out" molfoundry-observables eight-sites --out /out/eight-sites
```

For editable notebooks, mount a working copy of the checkout:

```sh
docker run --rm -p 127.0.0.1:8888:8888 -v "$PWD:/work" -w /work -e PYTHONPATH=/work/src:/work \
  --entrypoint jupyter molfoundry-observables lab \
  --ip=0.0.0.0 --no-browser --allow-root notebooks/paper_examples.ipynb
```

Open the localhost token URL printed by Jupyter. Docker notebooks are tested
on native Linux/amd64. On the tested macOS/ARM emulation, Jupyter kernels stall;
use the local installation above for notebooks. Emulated CLI durations are not
the paper's timings.

## Check and reproduce

```sh
python -m pytest tests notebooks
python scripts/reproduction/reproduce_paper.py run --profile smoke --output /tmp/mfobs-smoke
python scripts/reproduction/reproduce_paper.py run --profile paper --output /tmp/mfobs-paper
```

The paper profile computes current analyses and independent checks. It does not
rerun the retained FRAGGER measurements. [Claims and checks](docs/paper-claims.md)
identifies the scope; [benchmark instructions](scripts/comparison/ab_scaling/README.md)
explain external FRAGGER setup. No engine/compiler binaries are stored here.

The website and hosted notebook service are maintained separately. This repository
contains the software and reproducible notebooks.

## Repository layout

| Directory | Purpose |
|---|---|
| `src/` | Analysis library, shared algorithms and CLI support |
| `examples/` | Model inputs and example-specific reference data |
| `tutorials/` | Short readable analyses and their command-line entry point |
| `notebooks/` | Executable Jupyter notebooks and presentation helpers |
| `scripts/` | Reproduction, benchmarks and maintenance commands |
| `tests/` | Automated checks; recorded measurements live in `tests/data/` |
| `schema/` | JSON formats for model and chemical inputs |
| `docs/` | Setup, extension and validation guides |
| `.github/workflows/` | Automated checks on GitHub |

## License

Except where otherwise noted, the original software, notebooks and accompanying
documentation are licensed under [Apache 2.0](LICENSE). See [NOTICE.md](NOTICE.md).
Third-party components retain their own terms.
