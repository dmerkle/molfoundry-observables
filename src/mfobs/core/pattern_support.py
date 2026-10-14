"""Query-driven discovery of generators of a polynomial observable algebra.

The control algorithm contains no graph labels, motif lengths, candidate
patterns or target basis. A sound backend reduces actual projected generator
images and activity identities to polynomials in represented pattern counts.
Only factors occurring after signed cancellation are added to the queue.
Polynomial support closure is NOT linear expectation closure or jump closure.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Mapping, Protocol, Sequence

import sympy as sp


class UnsupportedReduction(ValueError):
    """A backend cannot justify this operation on its declared language."""


@dataclass
class Pattern:
    name: str
    rule: object
    kind: str
    metadata: dict = field(default_factory=dict)


@dataclass
class Image:
    expression: sp.Expr
    factors: Mapping[str, Pattern]
    evidence: dict = field(default_factory=dict)


class CountingBackend(Protocol):
    """Exact activities and generator images on one declared state family."""

    def activity(self, source) -> Image:
        ...

    def image(self, pattern: Pattern, source) -> Image:
        ...


@dataclass
class Support:
    status: str
    patterns: dict[str, Pattern]
    images: dict[str, dict[str, sp.Expr]]
    activities: dict[str, sp.Expr]
    origins: dict[str, dict]
    trace: list[dict]
    pending: list[str]
    reason: str | None = None

    def record(self) -> dict:
        """Serialize discovered counts and their exact generator identities."""
        return {
            "status": self.status,
            "selected_queries": [n for n, origin in self.origins.items()
                                 if origin["kind"] == "selected-query"],
            "patterns": {n: {"kind": p.kind, **p.metadata,
                              "origin": self.origins[n], "gml": p.rule.getGMLString()}
                         for n, p in self.patterns.items()},
            "activities": {n: polynomial_record(v) for n, v in self.activities.items()},
            "generator_images": {n: {r: polynomial_record(v) for r, v in row.items()}
                                 for n, row in self.images.items()},
            "origins": self.origins, "trace": self.trace,
            "pending": self.pending, "reason": self.reason,
        }


def polynomial_record(expression: sp.Expr) -> dict:
    expression = sp.expand(expression)
    names = sorted(str(s) for s in expression.free_symbols)
    symbols = [sp.Symbol(n) for n in names]
    if symbols:
        poly = sp.Poly(expression, *symbols, domain=sp.QQ)
        terms = [{"coefficient": str(c),
                  "powers": {n: k for n, k in zip(names, powers) if k}}
                 for powers, c in poly.terms() if c]
    else:
        if not expression.is_Rational:
            raise sp.PolynomialError("non-rational constant")
        terms = [] if expression == 0 else [{"coefficient": str(expression), "powers": {}}]
    return {"text": str(expression), "terms": terms}


def reduce_signed_rules(coefficients, representatives, reduce) -> Image:
    """Reduce a signed rule expression and collect its exact count factors."""
    expression, factors, terms = sp.Integer(0), {}, []
    for key, coefficient in sorted(coefficients.items()):
        if not coefficient:
            continue
        rule = representatives[key]
        try:
            image = reduce(rule)
        except UnsupportedReduction as exc:
            exc.term = {"key": key, "gml": rule.getGMLString(), "coefficient": coefficient}
            raise
        expression += coefficient * image.expression
        factors.update(image.factors)
        terms.append({"key": key, "gml": rule.getGMLString(), "coefficient": coefficient,
                      "polynomial": str(image.expression), "reduction": image.evidence})
    return Image(sp.expand(expression), factors, {"terms": terms})


def discover_support(
    rules: Sequence[object], seeds: Mapping[str, Pattern], *,
    counter: CountingBackend,
    maximum_patterns: int = 64, maximum_expansions: int = 64,
) -> Support:
    """Compute an exact polynomial dependency fixed point, or retain UNKNOWN.

    Rules have an ``identifier`` attribute. Budgets are checked before adding
    a new batch of patterns; the current obligation is preserved on stopping.
    Backend exceptions have a separate unsupported-reduction outcome.
    """
    if not seeds or maximum_patterns < len(seeds) or maximum_expansions < 0:
        raise ValueError("nonempty seeds and adequate nonnegative budgets required")
    identifiers = [r.identifier for r in rules]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("rule identifiers must be unique")
    patterns = dict(seeds)
    queue = deque(patterns)
    origins = {n: {"kind": "selected-query"} for n in patterns}
    images: dict[str, dict[str, sp.Expr]] = {}
    activities: dict[str, sp.Expr] = {}
    trace: list[dict] = []

    def finish(status: str, reason=None, current=None) -> Support:
        pending = ([current] if current and current not in queue else []) + list(queue)
        return Support(status, patterns, images, activities, origins, trace, pending, reason)

    def retain(image: Image, origin: dict) -> None:
        expr = sp.expand(image.expression)
        polynomial_record(expr)  # reject rational/non-polynomial results
        needed = sorted(str(s) for s in expr.free_symbols if str(s) not in patterns)
        missing = set(needed) - set(image.factors)
        if missing:
            raise UnsupportedReduction(f"unrepresented pattern factors: {sorted(missing)}")
        if len(patterns) + len(needed) > maximum_patterns:
            trace.append({**origin, "polynomial": polynomial_record(expr),
                          "pending_factors": needed, **image.evidence})
            raise OverflowError("pattern budget exhausted")
        for name in needed:
            pattern = image.factors[name]
            if pattern.name != name:
                raise UnsupportedReduction("factor name disagrees with its represented pattern")
            patterns[name] = pattern
            origins[name] = {**origin, "trace_index": len(trace)}
            queue.append(name)
        trace.append({**origin, "polynomial": polynomial_record(expr),
                      "new_patterns": needed, **image.evidence})

    current = None
    try:
        # Activities are derived from the supplied rules, not supplied counts.
        for r in rules:
            current = "activity:" + r.identifier
            result = counter.activity(r)
            retain(result, {"kind": "rule-activity", "rule": r.identifier})
            activities[r.identifier] = sp.expand(result.expression)
        expanded = 0
        while queue:
            if expanded >= maximum_expansions:
                return finish("UNKNOWN_EXPANSION_LIMIT", "expansion budget exhausted")
            current = queue.popleft()
            images[current] = {}
            for r in rules:
                result = counter.image(patterns[current], r)
                retain(result, {"kind": "generator-image", "pattern": current,
                                "rule": r.identifier})
                images[current][r.identifier] = sp.expand(result.expression)
            expanded += 1
        return finish("POLYNOMIAL_SUPPORT_CLOSED")
    except OverflowError as exc:
        return finish("UNKNOWN_PATTERN_LIMIT", str(exc), current)
    except (UnsupportedReduction, sp.PolynomialError) as exc:
        return finish("UNKNOWN_UNSUPPORTED_REDUCTION", str(exc), current)
