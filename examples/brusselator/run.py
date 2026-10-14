"""Brusselator–ROP through the public analysis API."""
from mfobs import BrusselatorROP, analyze
from mfobs.models.brusselator_pjc import selected_query_names


def run(*, linear_limit: int = 24, queries=None) -> dict:
    backend = BrusselatorROP(None if queries is None else selected_query_names(queries))
    return analyze(backend.specification(), backend=backend, linear_limit=linear_limit).record()
