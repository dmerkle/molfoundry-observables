"""Chlorination through the public analysis API."""
from mfobs import Chlorination, analyze


def run(variant="guarded", *, queries=None) -> dict:
    if queries is not None:
        raise ValueError("the neighbour-NAC model selects Z; custom queries are unsupported")
    backend = Chlorination(variant)
    return analyze(backend.specification(), backend=backend).record()
