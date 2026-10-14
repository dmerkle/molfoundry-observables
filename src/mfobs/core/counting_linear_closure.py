"""Algorithm 1 from conditional compositions and exact count reduction.

Keep graph expressions for subsequent compositions and reduce their counts
for exact span decisions. Neither support discovery nor event extraction is
used. The worklist itself is ``linear_closure.close_linear_span``.
"""
from collections import defaultdict
from fractions import Fraction

import sympy as sp

from mfobs.core.linear_closure import close_linear_span, UnsupportedImage
from mfobs.core.pattern_support import UnsupportedReduction, polynomial_record


def discover_linear_observables(rules, observables, *, counter,
                                maximum_dimension=24, input_definition_sha256=""):
    """Close the selected count functions under each rate-channel generator."""
    if len({rule.identifier for rule in rules}) != len(rules):
        raise ValueError("rule identifiers must be unique")
    channels = defaultdict(list)
    for rule in rules:
        channels[rule.rate].append(rule)
    graphs, expressions, forms = {}, {}, {}
    images = {}
    current = {}

    def vector(expression):
        return {tuple(sorted(term["powers"].items())): Fraction(term["coefficient"])
                for term in polynomial_record(expression)["terms"]}

    def polynomial(values):
        return sp.expand(sum(sp.Rational(c.numerator, c.denominator) *
                             sp.prod(sp.Symbol(n)**p for n, p in powers)
                             for powers, c in values.items()))

    def encode(values):
        return {"text": str(polynomial(values)),
                "degree": max((sum(p for _, p in powers) for powers in values), default=0), "terms": [
            {"monomial": dict(powers), "coefficient": str(c)}
            for powers, c in sorted(values.items()) if c]}

    def retain(rule, coefficient, target):
        reduced = counter.reduce(rule)
        if not reduced.expression:
            return
        # Equal reduced count functions have equal generator images on Ω.
        key = sp.srepr(sp.expand(reduced.expression))
        graphs.setdefault(key, rule)
        forms.setdefault(key, vector(reduced.expression))
        target[key] = target.get(key, Fraction(0)) + coefficient
        if not target[key]:
            del target[key]

    seeds = {}
    for name, pattern in observables.items():
        expressions[name] = {}
        terms = pattern.identity_terms() if hasattr(pattern, "identity_terms") else ((1, pattern.rule),)
        reduced = sp.Integer(0)
        for coefficient, rule in terms:
            retain(rule, Fraction(coefficient), expressions[name])
            reduced += coefficient * counter.reduce(rule).expression
        seeds[name] = vector(reduced)

    def image(source, channel):
        terms, evidence = {}, []
        try:
            for key, coefficient in expressions[source].items():
                for reaction in channels[channel]:
                    pair = (key, reaction.identifier)
                    if pair not in images:
                        calculated = counter.compose(graphs[key], reaction)
                        row = {}
                        for jump, value in calculated.observable_coefficients.items():
                            retain(calculated.observable_representatives[jump],
                                   Fraction(value) * reaction.per_match, row)
                        images[pair] = row, calculated.composition_orders
                    row, orders = images[pair]
                    for target, value in row.items():
                        terms[target] = terms.get(target, Fraction(0)) + coefficient * value
                    evidence.append({"source_rule": reaction.identifier, "composition_orders": orders})
        except UnsupportedReduction as exc:
            raise UnsupportedImage(str(exc)) from exc
        current.clear()
        current.update({key: value for key, value in terms.items() if value})
        total = {}
        for key, coefficient in terms.items():
            for monomial, value in forms[key].items():
                total[monomial] = total.get(monomial, Fraction(0)) + coefficient * value
        total = {key: value for key, value in total.items() if value}
        return total, {"generator_image": encode(total), "compositions": evidence}

    def designate(values, index):
        name = str(polynomial(values))
        return name if name not in expressions else f"O{index}"

    def added(name, values, source, channel, before, scale):
        terms = dict(current)
        for existing, coefficient in before.items():
            for key, value in expressions[existing].items():
                terms[key] = terms.get(key, Fraction(0)) - coefficient * value
        expressions[name] = {key: value / scale for key, value in terms.items() if value}

    result = close_linear_span(seeds=seeds, channels=channels, image=image,
        encode=encode, designate=designate, on_add=added,
        maximum_dimension=maximum_dimension)
    return {"schema_version": 2, "status": result.status,
            "analysis": "linear-observable-discovery",
            "input_definition_sha256": input_definition_sha256,
            "basis": list(result.basis),
            "basis_expressions": {name: encode(values) for name, values in result.basis.items()},
            "origins": result.origins, "steps": result.steps,
            "pending": result.pending, "obligations": result.obligations,
            "dimension": len(result.basis), "maximum_dimension": maximum_dimension,
            "scope": "exact conditional compositions and count reductions on Ω"}
