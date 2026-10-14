"""Rémy trees through the public analysis API."""
from mfobs import RemyTrees, analyze
from mfobs.models.remy.inputs import query_order


def run(*, linear_limit: int = 16, queries=None) -> dict:
    backend = RemyTrees(3 if queries is None else query_order(queries))
    return analyze(backend.specification(), backend=backend, linear_limit=linear_limit).record()
