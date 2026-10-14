"""Generic no-reduction adapter for new small rule sets.

Canonical equality and direct edge/NoEdge contradiction elimination are the
only reductions. This may not terminate on an infinite
hierarchy: the shared exploration budget returns UNKNOWN, not a false proof.
No state-language invariant or chemistry knowledge is assumed by this adapter.
"""
from fractions import Fraction
from mfobs.core.observables import parse_observable
from mfobs.core.conditions import (
    adjacency_contradictions,
    engine_conditions,
    identity_rule,
)


def noedge_contradictions(rule):
    """Pairs simultaneously required to be edges and not edges on the input.

    This is an exact zero test on the supported undirected simple-graph
    semantics, independent of any supplied state language. A negative result
    is NOT a satisfiability certificate: no general condition solver is used.
    """
    parsed = parse_observable(rule)
    positive_pairs = {(u, v) for u, v, _ in parsed.edges}
    return tuple(sorted(positive_pairs & set(parsed.no_edges)))


def identity_contradictions(rule):
    """Certified direct contradictions for one supported identity pattern."""

    reasons = []
    pairs = noedge_contradictions(rule)
    if pairs:
        reasons.append({
            'kind': 'positive-edge-and-NoEdge',
            'pairs': [list(pair) for pair in pairs],
        })
    for contradiction in adjacency_contradictions(rule):
        reasons.append({
            'kind': 'incompatible-adjacency-conditions',
            **contradiction,
        })
    return tuple(reasons)


class DanglingNotCertified(ValueError):
    """A DPO jump closure deletes a vertex whose dangling condition is not certified."""


def jump_input(rule):
    """Parsed input ``I`` of a DPO jump closure ``I <- K -> I`` and its vertices ``I \\ K``.

    ``dpo_jump_closure`` stores a deleted vertex on the left and recreates it
    under a fresh id on the right.  Raises ValueError unless the right side is
    the left side up to that renaming, i.e. unless ``rule`` is identity-like.
    """

    from collections import Counter
    from mfobs.core.rules import restores_input, rule_vertex_ids

    deleted, created = rule_vertex_ids(rule, 'deleted'), rule_vertex_ids(rule, 'created')
    left = {v.id: v.leftLabel for v in rule.vertices}
    right = {v.id: v.rightLabel for v in rule.vertices}
    if (Counter(left[v] for v in deleted)
            != Counter(right[v] for v in created)):
        raise ValueError('not a DPO jump closure: deleted and recreated vertices differ')
    if not restores_input(rule):
        raise ValueError('not a DPO jump closure: right side is not a recreated left side')
    return parse_observable(rule), deleted


def certified_input_identity(rule, certify=None):
    """Reduce a DPO jump closure to the identity counting its input ``I``, if sound.

    Without deleted vertices the closure counts matches of ``I``.  A deleted
    vertex adds its DPO dangling condition (no host edge outside the match).
    ``certify(graph, deleted)`` must show that this holds on every state of the
    declared language, raising :class:`DanglingNotCertified` otherwise.
    """

    parsed, deleted = jump_input(rule)
    if deleted:
        if certify is None:
            raise DanglingNotCertified(
                f'deleted vertices {list(deleted)} need a dangling certificate '
                'for the declared state language')
        from mfobs.core.graph import make_graph
        certify(make_graph(parsed.labels, parsed.edges), deleted)
    return _identity_of(rule, parsed)


def normalized_identity_rule(rule):
    """Rebuild one supported identity with duplicate conditions removed."""

    parsed = parse_observable(rule)
    if not parsed.jump_identity:
        raise ValueError('identity normalization supports only direct conditions')
    return _identity_of(rule, parsed)


def _identity_of(rule, parsed):
    from mfobs.core.graph import make_graph

    if parsed.generic_conditions not in {(), ('matchConstraintsAdjacency',)}:
        raise ValueError('identity normalization supports only direct conditions')
    graph = make_graph(parsed.labels, parsed.edges)
    adjacency = tuple(sorted(set(engine_conditions(rule))))
    return identity_rule(graph, 'normalized_identity', tuple(sorted(set(parsed.no_edges))), adjacency)


from mfobs.core.synthesis import DynamicBasisAdapter, Reduction
from mfobs.core.linear_closure import RationalSpan


class Adapter(DynamicBasisAdapter):
    semantic_contract_required = False

    def __init__(self, inputs):
        self.inputs = inputs
        self.span = RationalSpan()
        self.forms, self.origins, self.jump_terms, self.zero_jump_terms = {}, {}, {}, {}

    def reduce_pattern(self, rule, key):
        contradictions = noedge_contradictions(rule)
        if contradictions:
            self.zero_jump_terms.setdefault(key, {
                'reason': 'positive-edge-and-NoEdge-on-same-input-pair',
                'contradictory_pairs': [list(pair) for pair in contradictions],
                'gml': rule.getGMLString()})
            return Reduction('zero', {})
        return Reduction('exact', {key: Fraction(1)})



    def encode_form(self, form):
        return {key: str(c) for key, c in form.items()}

    def finalize(self, state, composer, complete):
        return {'status': 'CLOSED_CANONICAL_JUMP_SPANS',
                'proof_status': 'Exact channel residual closure under canonical identity equality and direct edge/NoEdge contradiction elimination.',
                'jump_terms': self.jump_terms, 'zero_jump_terms': self.zero_jump_terms,
                'initial_values': {}}
