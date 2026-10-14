"""Public MolFoundry composition with state-language constraints and predicates.

MolFoundry enumerates overlaps, checks DPO admissibility and transports
conditions. Local structural theories are passed as MolFoundry incidence
constraints, which the engine checks on all three graphs of an overlap and
uses to skip partial overlaps no extension can satisfy. Python predicates
check only what is not local to one vertex. Every surviving overlap is
grouped; canonically equal products retain their multiplicity.
"""
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field

import molfoundry as mf


class CompositionLimit(RuntimeError):
    """No partial composition counter may be used after this exception."""


@dataclass
class CompositionDiagnostics:
    # The public iterator does not expose internal gluing/condition rejections.
    candidate_overlaps: int | None = None
    null_compositions: int | None = None
    condition_rejected: int | None = None
    forbidden_overlap_rejected: int = 0
    forbidden_input_rejected: int = 0
    forbidden_output_rejected: int = 0
    accepted_raw: int = 0
    forbidden_reasons: dict = field(default_factory=dict)
    predicate_checks: dict = field(default_factory=lambda: dict(overlap=0, input=0, output=0))

    def record(self):
        return asdict(self)


class Representatives(Mapping):
    """Keep one GML representative per canonical class.

    Rebuild a rule on lookup instead of retaining every engine rule in memory.
    """

    def __init__(self, gml=None):
        self._gml = dict(gml or {})

    def add(self, key, rule):
        if key not in self._gml:
            self._gml[key] = rule.getGMLString()

    def gml(self, key):
        return self._gml[key]

    def __getitem__(self, key):
        return mf.Rule.fromGMLString(self._gml[key], add=False)

    def __iter__(self):
        return iter(self._gml)

    def __len__(self):
        return len(self._gml)

    def __contains__(self, key):
        return key in self._gml

    def __or__(self, other):
        """Both sides' classes, without rebuilding any rule."""
        return Representatives({**self._gml, **other._gml})


def rule_labels(*rules):
    """Vertex and edge labels on either side of the rules; every graph of their
    compositions is labelled from these finite sets."""
    vertices, edges = set(), set()
    for rule in rules:
        for v in rule.vertices:
            vertices.update(label for label in (v.leftLabel, v.rightLabel) if label is not None)
        for e in rule.edges:
            edges.update(label for label in (e.leftLabel, e.rightLabel) if label is not None)
    return vertices, edges


def group_compositions(first, second, *, violation=None, constraints=None, connected=False,
                       maximum=False, include_empty=True, maximum_intermediate_graphs=None,
                       rule_rejection=None, on_rule=None):
    """Return diagnostics, condition-sensitive canonical multiplicities and
    one representative per class (see ``Representatives``).

    MolFoundry checks ``constraints`` before ``violation`` sees a graph. The
    optional limit counts intermediate graphs that pass the constraints, not
    MolFoundry's internal search steps. It cannot bound work inside next().
    """
    diagnostics = CompositionDiagnostics()

    def predicate(phase):
        def accepts(graph):
            diagnostics.predicate_checks[phase] += 1
            if (phase == 'overlap' and maximum_intermediate_graphs is not None
                    and diagnostics.predicate_checks[phase] > maximum_intermediate_graphs):
                raise CompositionLimit('COMPOSITION_CANDIDATE_LIMIT: intermediate graph limit exceeded; no partial count')
            if violation is None:
                return True
            reason = violation(graph)
            if reason is None:
                return True
            attribute = {'overlap': 'forbidden_overlap_rejected',
                         'input': 'forbidden_input_rejected',
                         'output': 'forbidden_output_rejected'}[phase]
            setattr(diagnostics, attribute, getattr(diagnostics, attribute) + 1)
            key = f'{phase}:{reason}'
            diagnostics.forbidden_reasons[key] = diagnostics.forbidden_reasons.get(key, 0) + 1
            return False
        return accepts

    if violation is not None:
        predicates = dict(intermediatePredicate=predicate('overlap'),
                          leftPredicate=predicate('input'), rightPredicate=predicate('output'))
    elif maximum_intermediate_graphs is not None:
        predicates = dict(intermediatePredicate=predicate('overlap'))
    else:
        predicates = {}
    evaluator = mf.RCEvaluator(constraints=constraints, **predicates)
    expression = first * mf.rcCommon(connected=connected, maximum=maximum,
                                     includeEmpty=include_empty) * second
    counts, representatives = Counter(), Representatives()
    iterator = evaluator.iterate(expression, onlyUnique=False)
    try:
        for rule in iterator:
            try:
                reason = rule_rejection(rule) if rule_rejection else None
                if on_rule is not None:
                    on_rule(rule, reason)
                if reason:
                    diagnostics.forbidden_input_rejected += 1
                    diagnostics.forbidden_reasons[reason] = diagnostics.forbidden_reasons.get(reason, 0) + 1
                else:
                    key = rule.canonKey
                    counts[key] += 1
                    representatives.add(key, rule)
            finally:
                evaluator.createdRules.clear()
    finally:
        iterator.close()
    diagnostics.accepted_raw = sum(counts.values())
    return diagnostics, counts, representatives


def witnessed_compositions(first, second, *, include_empty=True):
    """Yield each public composite and its first-right to second-left overlap."""
    overlap = ()
    def remember(graph):
        nonlocal overlap
        overlap = tuple((v.firstId, v.secondId) for v in graph.vertices
                        if v.firstId is not None and v.secondId is not None)
        return True
    evaluator = mf.RCEvaluator(intermediatePredicate=remember)
    expression = first * mf.rcCommon(connected=False, includeEmpty=include_empty) * second
    iterator = evaluator.iterate(expression, onlyUnique=False)
    try:
        for rule in iterator:
            yield rule, overlap
            evaluator.createdRules.clear()
    finally:
        iterator.close()
