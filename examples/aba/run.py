"""Export the growing-path analysis for a literal ABA query."""
import json
from pathlib import Path
from typing import Mapping

from mfobs import analyze, ABA
from mfobs.models.aba import counting as terminal


QUERY = Path(__file__).with_name("query.json")


def run(query: Mapping | None = None, *, diagnostic_timings: bool = False,
            queries=None) -> dict:
    query = dict(query) if query is not None else json.loads(QUERY.read_text(encoding="utf-8"))
    motifs = terminal.query_motifs(queries) if queries is not None else (query["motif"],)
    backend = ABA(motifs, alphabet=tuple(query["alphabet"]),
                  gated=bool(query.get("gate", query.get("gated", False))))
    return analyze(
        backend.specification(), backend=backend,
        linear_limit=int(query["linear_limit"]),
        support_pattern_limit=int(query["support_pattern_limit"]),
        support_expansion_limit=int(query["support_expansion_limit"]),
        diagnostic_timings=diagnostic_timings,
    ).record()


