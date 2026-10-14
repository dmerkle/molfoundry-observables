"""Independent host matching for exported linear combinations of patterns."""
from fractions import Fraction

from mfobs.core.conditions import parse_observable
from mfobs.core.engine import mf

from mfobs.core.mixtures import fast_count


def parse_expressions(result):
    expressions = result.get('observable_expressions', {
        name: [{'coefficient': '1', 'gml': item['gml']}]
        for name, item in result['identity_patterns'].items()})
    if set(expressions) != set(result['basis']):
        raise ValueError('Every basis coordinate needs an explicit observable expression')
    parsed = {}
    answer = {}
    for name, terms in expressions.items():
        answer[name] = []
        for term in terms:
            gml = term['gml']
            if gml not in parsed:
                rule = mf.Rule.fromGMLString(gml, add=False)
                pattern, _ = parse_observable(rule)
                if not pattern.jump_identity or pattern.generic_conditions:
                    raise ValueError('Host checks require supported conditional identity patterns')
                parsed[gml] = rule
            answer[name].append((Fraction(term['coefficient']), parsed[gml]))
    return answer


def count(expression, host):
    return sum(c * fast_count(p, host)
               for c, p in expression)
