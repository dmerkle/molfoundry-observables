"""Ring8 through the public analysis API."""
from mfobs import Ring8, analyze
from mfobs.inputs import with_observable_queries


def run(*, worklist_limit: int = 1000, queries=None, include_evidence=False) -> dict:
    backend = Ring8()
    if queries is not None:
        backend = Ring8(manifest=with_observable_queries(backend.inputs(), queries))
    return analyze(backend.specification(), backend=backend, worklist_limit=worklist_limit,
                   include_evidence=include_evidence).record()


def run_with_evidence(*, worklist_limit: int = 1000, queries=None) -> tuple[dict, dict]:
    result = run(worklist_limit=worklist_limit, queries=queries, include_evidence=True)
    raw = result.pop("linear_evidence")
    return result, raw
