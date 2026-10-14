"""Fresh shared analyses for the CLI and numerical regression checks."""
import gzip
import json
from pathlib import Path
import pytest


@pytest.fixture(scope="session")
def aba_records():
    from mfobs import ABA, analyze
    return {name: analyze((backend := ABA.motif("ABA", gated=gated)).specification(),
                          backend=backend, linear_limit=24).record()
            for name, gated in (("aba", False), ("aba-gated", True))}


@pytest.fixture(scope="session")
def chlorination_records():
    from mfobs import Chlorination, analyze
    return {name: analyze((backend := Chlorination(name)).specification(), backend=backend).record()
            for name in ("guarded", "ablated")}


@pytest.fixture(scope="session")
def remy_record():
    from mfobs import RemyTrees, analyze
    backend = RemyTrees()
    return analyze(backend.specification(), backend=backend).record()


def regression_record(name):
    with gzip.open(Path(__file__).parent / "tests/data" / (name + ".json.gz"), "rt", encoding="utf-8") as stream:
        return json.load(stream)


@pytest.fixture(scope="session")
def ring_composition():
    from mfobs.core.worklist import run
    from mfobs.inputs import load
    return run(load(Path(__file__).parent / "examples/ring8/model.json"))


@pytest.fixture(scope="session")
def eight_composition():
    from mfobs.core.worklist import run
    from mfobs.inputs import load
    return run(load(Path(__file__).parent / "examples/eight-slots/model.json"))
