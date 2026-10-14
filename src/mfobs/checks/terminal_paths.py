"""Independent word oracle for the terminal-path and molecular ROP analyses."""
from collections import Counter
from itertools import product
from typing import Mapping

import sympy as sp

from mfobs.core.pattern_support import Support, polynomial_record
from mfobs.models.aba.events import compile_uniform_events

def coordinate_values(patterns, word: str, pools: Mapping[str, int], opened=True) -> dict[str, int]:
    """Independent word evaluation, used only after discovery for validation."""
    values = {}
    for name, pattern in patterns.items():
        w = pattern.metadata.get("word", "")
        if pattern.kind == "pool":
            value = pools[w]
        elif pattern.kind == "total":
            value = pools[w] + word.count(w)
        elif pattern.kind == "motif":
            value = sum(word[i:i+len(w)] == w for i in range(len(word)-len(w)+1))
        elif pattern.kind == "suffix":
            value = int(word.endswith(w))
        elif pattern.kind == "gate":
            value = int(opened)
        else:
            raise ValueError("unknown count meaning")
        values[name] = value
    return values


def validate_words(support: Support, sources, *, maximum_length=6, jump_model=None, event_compiler=None) -> dict:
    """Compare complete jump/rate distributions, not merely their first images."""
    from fractions import Fraction
    event_result = (event_compiler or compile_uniform_events)(support, sources)
    if event_result["status"] != "FINITE_EVENTS":
        return {"status": "NOT_EXECUTED", "reason": event_result["status"]}
    events = event_result["events"]
    names = list(support.patterns)
    alphabet = tuple(r.symbol for r in sources if r.kind == "feed")
    gated = any(r.kind == "bind" for r in sources)

    def evaluator(polynomial):
        terms = [(Fraction(t["coefficient"]), tuple(t["powers"].items()))
                 for t in polynomial["terms"]]
        return lambda v: sum(c * sp.prod(v[n] ** p for n, p in powers) for c, powers in terms)

    image_evals = {n: {r.identifier: evaluator(polynomial_record(support.images[n][r.identifier]))
                       for r in sources} for n in names}
    activity_evals = {r.identifier: evaluator(polynomial_record(support.activities[r.identifier])) for r in sources}
    event_evals = [(e, evaluator(e["activity"])) for e in events]
    states, scalar_identities, distributions = 0, 0, 0
    for length in range(maximum_length + 1):
        for letters in product(alphabet, repeat=length):
            word = "".join(letters)
            for pool_vector in product((0, 1, 2), repeat=len(alphabet)):
                pools = dict(zip(alphabet, pool_vector))
                for opened in ((False, True) if gated else (True,)):
                    before = coordinate_values(support.patterns, word, pools, opened)
                    states += 1
                    for r in sources:
                        target_word, target_pool, target_open = word, dict(pools), opened
                        if r.kind == "feed":
                            multiplicity = 1
                            target_pool[r.symbol] += 1
                        elif r.kind == "append":
                            multiplicity = pools[r.symbol] * int(opened)
                            if multiplicity:
                                target_word += r.symbol
                                target_pool[r.symbol] -= 1
                        elif r.kind == "bind":
                            multiplicity, target_open = int(opened), False
                        else:
                            multiplicity, target_open = int(not opened), True
                        after = coordinate_values(support.patterns, target_word, target_pool, target_open)
                        delta = tuple(after[n] - before[n] for n in names)
                        assert activity_evals[r.identifier](before) == multiplicity
                        for index, n in enumerate(names):
                            assert image_evals[n][r.identifier](before) == multiplicity * delta[index], (word, n, r.identifier)
                            scalar_identities += 1
                        observed = Counter()
                        for event, ev in event_evals:
                            if event["source_rule"] != r.identifier:
                                continue
                            weight = ev(before)
                            assert weight >= 0
                            if weight:
                                observed[tuple(event["jump"].get(n, 0) for n in names)] += weight
                        expected = Counter({delta: multiplicity}) if multiplicity and any(delta) else Counter()
                        assert observed == expected, (word, pools, opened, r.identifier, observed, expected)
                        if jump_model is not None:
                            converted = Counter()
                            for event in jump_model.events:
                                if event.rule != r.identifier:
                                    continue
                                weight = jump_model.ring.evaluate(event.propensity, before)
                                assert weight >= 0
                                if weight:
                                    converted[event.delta] += weight
                            assert tuple(jump_model.ring.variables) == tuple(names)
                            assert converted == expected, (word, pools, opened, r.identifier, converted, expected)
                        distributions += 1
    return {"status": "PASS", "maximum_word_length": maximum_length, "pool_values": [0, 1, 2],
            "states": states, "scalar_generator_identities": scalar_identities,
            "channel_jump_distributions": distributions,
            "adapted_model_jump_distributions": distributions if jump_model is not None else 0,
            "scope": "bounded independent word oracle after all symbolic discovery and event extraction"}


def result_components(record, *, molecular=False):
    """Read computed patterns and events for an independent regression check."""
    from dataclasses import replace
    from types import SimpleNamespace
    import molfoundry as mf
    from mfobs.core.pattern_support import Pattern
    from mfobs.core.polynomial_jump_closure import JumpEvent, PolynomialRing
    from mfobs.models.aba.inputs import source_rules
    from mfobs.models.molecular_rop import inputs as rop

    data = record["configuration"]["model"]
    if molecular:
        selected = rop.RopConfiguration(data["identifier"], data["title"],
            tuple(data["alphabet"]), data["gate"], tuple(tuple(w) for w in data["motifs"]), data["purpose"])
    else:
        selected = SimpleNamespace(identifier=data["identifier"], title=data["title"],
            alphabet=tuple(data["alphabet"]), gated=data["gate"], motifs=tuple(tuple(w) for w in data["motifs"]))
    raw = record["support_discovery"]
    patterns = {n: Pattern(n, mf.Rule.fromGMLString(p["gml"], add=False), p["kind"],
                          {k: v for k, v in p.items() if k not in ("kind", "gml", "origin")})
                for n, p in raw["patterns"].items()}
    symbols = {n: sp.Symbol(n) for n in patterns}
    read = lambda p: sp.sympify(p["text"], locals=symbols)
    support = Support(raw["status"], patterns,
        {n: {r: read(p) for r, p in row.items()} for n, row in raw["generator_images"].items()},
        {r: read(p) for r, p in raw["activities"].items()},
        raw["origins"], raw["trace"], raw["pending"], raw["reason"])
    constructors = (rop.source_rules(selected) if molecular else
                    source_rules(selected.alphabet, selected.gated))
    sources = tuple(replace(r, rule=mf.Rule.fromGMLString(raw["source_rules"][r.identifier]["gml"], add=False),
                            rate=raw["source_rules"][r.identifier]["rate"]) for r in constructors)
    model = None
    if raw["events"]["status"] == "FINITE_EVENTS":
        ring = PolynomialRing(tuple(patterns))
        model = SimpleNamespace(ring=ring, events=tuple(
            JumpEvent(e["identifier"], e["rate"], tuple(e["jump"].get(n, 0) for n in patterns),
                      ring.polynomial({tuple(t["monomial"].get(n, 0) for n in patterns): t["coefficient"]
                                       for t in e["activity"]["terms"]}), e["source_rule"], e["predicate"])
            for e in record["composition"]["event_classes"]))
    return selected, support, sources, model


def validate_result(record, *, maximum_length=6):
    """Compare the public ABA result with concrete word transitions."""
    _, support, sources, model = result_components(record)
    return validate_words(support, sources, maximum_length=maximum_length, jump_model=model)
