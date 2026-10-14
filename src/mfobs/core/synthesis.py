"""Graph composition and exact counting for Algorithm 1.

Backends supply reductions and preservation evidence. The only closure loop is
``linear_closure.close_linear_span``. Here generator images and residual graph
expressions are assembled without choosing a stronger pattern/monomial basis.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Iterable, Literal

from mfobs.core.exact_coefficients import add_vectors
from mfobs.core.observables import parse_observable
from mfobs.inputs import InputError


@dataclass(frozen=True)
class Reduction:
    """Exact semantic reduction of one identity/jump observable."""

    status: Literal['exact', 'zero', 'unknown']
    form: Any
    witness: Any = None
    reason: str = ''
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.status not in {'exact', 'zero', 'unknown'}:
            raise InputError(f'Invalid reduction status: {self.status}')
        if self.status == 'unknown' and not self.reason:
            raise InputError('An UNKNOWN reduction needs an explicit reason')


@dataclass(frozen=True)
class Candidate:
    """Inspectable pattern that may represent a new basis coordinate."""

    key: str
    rule: object
    form: Any
    origin: dict
    witness: Any = None
    metadata: dict = field(default_factory=dict)


class DynamicBasisAdapter:
    """Template method implementing the shared synthesis algorithm.

    Subclasses implement only the state-language semantics: reduction, form serialization, candidate generation, and certificates.
    MolFoundry composition, the worklist, naming, rate separation, and dynamic
    span extension live here and in :mod:`mfobs.core.worklist`.
    """

    semantic_contract_required = True

    def __init__(self, inputs, span, semantic_contract):
        self.inputs = inputs
        self.span = span
        semantic_contract.validate_rules(
            rule['id'] for rule in inputs.data['rules']
        )
        self.semantic_contract = semantic_contract
        self.forms = {}
        self.origins = {}
        self.jump_terms = {}
        self.discovery_decisions = []
        self._candidates: dict[str, Candidate] = {}

    # ---- State-language semantics supplied by the application adapter. ----
    def reduce_pattern(self, rule, key: str) -> Reduction:  # pragma: no cover - abstract
        raise NotImplementedError

    def zero_form(self):
        return {}

    def add_form(self, total, form, coefficient):
        return add_vectors(total, form, coefficient)

    def encode_form(self, form):  # pragma: no cover - abstract
        raise NotImplementedError

    def form_support(self, form):
        return set(form)

    def meaning(self, form):
        return str(form)

    def seed_rule(self, supplied_rule, reduction: Reduction):
        return supplied_rule

    def seed_meaning(self, name: str, reduction: Reduction):
        return self.meaning(reduction.form)

    def new_observable_name(self, state):
        return f'O{len(state.patterns)}'

    def on_seed_added(self, name: str, reduction: Reduction):
        pass

    def on_candidate_added(self, name: str, candidate: Candidate):
        pass

    def candidate_patterns(
        self,
        rule,
        key: str,
        reduction: Reduction,
    ) -> Iterable[Candidate]:
        """Return actual jump patterns or justified projections thereof."""
        yield Candidate(
            key,
            rule,
            reduction.form,
            {'kind': 'actual-computed-jump', 'source_jump': key},
            reduction.witness,
        )

    def after_initialize(self):
        """Optionally register language-specific candidates such as identity(1)."""

    def seed_origin(self, name: str, reduction: Reduction) -> dict:
        return {'kind': 'input-seed'}

    def candidate_rank(self, candidate: Candidate, residual, remainder):
        direct = self.form_support(remainder) <= self.form_support(residual)
        parsed = parse_observable(candidate.rule)
        return (
            not direct,
            len(parsed.labels),
            len(parsed.no_edges) + len(parsed.generic_conditions),
            len(parsed.edges),
            candidate.key,
        )

    def format_coefficient(self, value):
        return str(value)

    # ---------------------- Shared executed algorithm. ----------------------
    def register_candidate(self, candidate: Candidate):
        is_new = candidate.key not in self._candidates
        is_nonzero = bool(self.form_support(candidate.form))
        if is_new and is_nonzero:
            self._candidates[candidate.key] = candidate

    def initialize(self, state, composer):
        self._candidates = {}
        self.discovery_decisions = []
        available = self.inputs.engine_observables()
        for name in self.inputs.data['initial_observables']:
            rule = available[name]
            reduced = self.reduce_pattern(rule, rule.canonKey)
            if reduced.status == 'unknown':
                raise InputError(
                    f'Initial observable {name} cannot be counted: '
                    f'{reduced.reason}'
                )
            if reduced.status == 'zero' or not self.form_support(reduced.form):
                raise InputError(
                    f'Initial observable {name} is zero on the declared language'
                )
            if not self.span.insert(name, reduced.form):
                raise InputError(
                    'Initial observables must be linearly independent on the '
                    'declared language'
                )
            self.forms[name] = reduced.form
            self.origins[name] = self.seed_origin(name, reduced)
            state.add(
                name,
                self.seed_rule(rule, reduced),
                self.seed_meaning(name, reduced),
            )
            self.on_seed_added(name, reduced)
        self.after_initialize()

    def run_closure(self, state, composer, *, max_expansions, maximum_dimension):
        from mfobs.core.linear_closure import close_linear_span, primitive, UnsupportedImage
        import sympy as sp

        self.discovery_decisions = []
        self.expressions = {
            name: {p.rule.canonKey: (Fraction(1), p.rule)}
            for name, p in state.patterns.items()
        }
        channels = {}
        for reaction in self.inputs.data['rules']:
            channels.setdefault(reaction['rate'], []).append(reaction)
        # Independent per-rule rows need not lie in the channel-closed span.
        # They remain in the trace; the exported ODE is explicitly by channel.
        state.equations = {}
        cache = {}
        images = {}
        forms_by_rule = {}
        current_image = {}

        def expression_add(left, right, factor=1):
            out = dict(left)
            for key, (value, rule) in right.items():
                coefficient = sp.cancel(out.get(key, (0, None))[0] + factor * value)
                if coefficient:
                    out[key] = (coefficient, rule)
                else:
                    out.pop(key, None)
            return out

        def image(source, channel):
            state.current_depth = state.patterns[source].depth
            total, terms = self.zero_form(), {}
            records = []
            for reaction in channels[channel]:
                per_rule = self.zero_form()
                for key, (coefficient, rule) in self.expressions[source].items():
                    cache_key = (key, reaction['id'])
                    if cache_key not in cache:
                        commutator, record = composer.pair(source, rule, reaction)
                        zero = self.pattern_derivative_zero_certificate(rule, reaction, commutator)
                        form, projected = self.zero_form(), {}
                        trace = []
                        if zero is not None:
                            record['derivative_zero_certificate'] = zero
                        else:
                            for jump_key, factor in sorted(commutator.observable_coefficients.items()):
                                jump = commutator.observable_representatives[jump_key]
                                reduced = self.reduce_pattern(jump, jump_key)
                                self.jump_terms.setdefault(jump_key, {
                                    'gml': jump.getGMLString(), 'reduction': reduced.witness})
                                trace.append({'key': jump_key, 'coefficient': str(factor)})
                                if reduced.status == 'unknown':
                                    raise UnsupportedImage(f"{reaction['id']}/{jump_key}: {reduced.reason}")
                                if reduced.status == 'zero':
                                    continue
                                form = self.add_form(form, reduced.form, factor)
                                representative = self.seed_rule(jump, reduced)
                                projected = expression_add(projected, {
                                    representative.canonKey: (factor, representative)})
                                # A collinear actual pattern is only a display/scale choice.
                                for candidate in self.candidate_patterns(jump, jump_key, reduced):
                                    self.register_candidate(candidate)
                        record.update(terms=trace, exact_derivative_form=self.encode_form(form))
                        cache[cache_key] = form, projected, record
                    form, projected, record = cache[cache_key]
                    factor = coefficient * Fraction(reaction['per_match'])
                    per_rule = self.add_form(per_rule, form, factor)
                    terms = expression_add(terms, projected, factor)
                    records.append({'composition_record': composer.records.index(record),
                                    'source_coefficient': str(coefficient),
                                    'per_match': reaction['per_match']})
                forms_by_rule[source, reaction['id']] = per_rule
                total = self.add_form(total, per_rule, 1)
            images[source, channel] = terms
            current_image.clear()
            current_image.update(total)
            return total, {'composition_terms': records}

        selected = {}
        self.span = type(self.span)()
        def normalize(residual):
            # Prefer a simple computed pattern only when its remainder spans
            # exactly the same one-dimensional extension (checked again by core).
            choices = []
            for candidate in self._candidates.values():
                remainder, _ = self.span.reduce(candidate.form)
                if set(remainder) != set(residual):
                    continue
                key = min(residual)
                scale = sp.cancel(residual[key] / remainder[key])
                if not self.add_form(residual, remainder, -scale):
                    choices.append((self.candidate_rank(candidate, residual, remainder), candidate, scale))
            if choices:
                _, candidate, scale = min(choices, key=lambda row: row[0])
                selected['candidate'] = candidate
                return candidate.form, scale
            selected.clear()
            if any(getattr(c, 'free_symbols', ()) for c in current_image.values()):
                # Over Q(parameters), the whole image has the same remainder.
                # Keeping it avoids introducing poles merely to make a pivot 1.
                selected['whole_image'] = True
                return current_image.copy(), Fraction(1)
            return primitive(residual)

        def simplify_expression(vector, fallback):
            # This auxiliary span only chooses a cheaper *representation* of
            # an already fixed observable. Its patterns never enter the FIFO
            # or the discovered basis. Every identity uses exact backend forms.
            options = []
            for index, candidate in enumerate(self._candidates.values()):
                rule = candidate.rule
                options.append((len(parse_observable(rule).labels), f'C{index}',
                                candidate.form, {rule.canonKey: (Fraction(1), rule)}))
            for name, form in self.forms.items():
                expression = self.expressions[name]
                cost = sum(len(parse_observable(rule).labels) for _, rule in expression.values())
                options.append((cost, f'B{name}', form, expression))
            from mfobs.core.linear_closure import RationalSpan
            def coefficients(form):
                # Require rational coefficients in this representation change.
                # Symbolic parameter monomials are formal components, so a
                # simplification cannot introduce a new parameter denominator.
                result = {}
                for key, value in form.items():
                    for monomial, coefficient in sp.expand(value).as_coefficients_dict().items():
                        if coefficient:
                            result[key, str(monomial)] = Fraction(coefficient)
                return result
            representation_span = RationalSpan()
            target_vector = coefficients(vector)
            representatives = {}
            for _, key, form, expression in sorted(options, key=lambda item: item[:2]):
                if not representation_span.insert(key, coefficients(form)):
                    continue
                representatives[key] = expression
                remainder, coordinates = representation_span.reduce(target_vector)
                if not remainder:
                    answer = {}
                    for target, coefficient in coordinates.items():
                        answer = expression_add(answer, representatives[target], coefficient)
                    return answer
            return fallback

        def added(name, vector, source, channel, before, scale):
            expr = images[source, channel]
            if not selected.get('whole_image'):
                for target, coefficient in before.items():
                    expr = expression_add(expr, self.expressions[target], -coefficient)
            expr = expression_add({}, expr, 1 / scale)
            candidate = selected.get('candidate')
            if candidate is not None:
                expr = {candidate.rule.canonKey: (Fraction(1), candidate.rule)}
                self.on_candidate_added(name, candidate)
            else:
                expr = simplify_expression(vector, expr)
            self.expressions[name] = expr
            self.forms[name] = vector
            self.origins[name] = {'kind': 'computed-channel-residual',
                                  'source': source, 'channel': channel}
            singleton = next(iter(expr.values())) if len(expr) == 1 else None
            rule = singleton[1] if singleton is not None and singleton[0] == 1 else None
            state.add(name, rule, self.meaning(vector))

        result = close_linear_span(
            seeds=self.forms.copy(), channels=channels, image=image,
            span=self.span, normalize=normalize, encode=self.encode_form,
            designate=lambda v, i: self.new_observable_name(state), on_add=added,
            on_expansion=state.record_expansion,
            maximum_dimension=maximum_dimension, maximum_expansions=max_expansions,
        )
        self.forms = result.basis
        self.discovery_decisions = result.steps
        state.queue.clear()
        state.queue.extend(result.pending)
        state.obligations = result.obligations
        state.channel_equations = {}
        for step in result.steps:
            if 'coordinates' in step:
                state.channel_equations.setdefault(step['source'], {})[step['channel']] = step['coordinates']
        # Retain per-rule coefficients only when they are actually representable.
        final_span = type(self.span)()
        for name, form in self.forms.items():
            final_span.insert(name, form)
        for (source, reaction), form in forms_by_rule.items():
            residual, coordinates = final_span.reduce(form)
            if not residual:
                state.equations.setdefault(source, {})[reaction] = {
                    k: self.format_coefficient(v) for k, v in coordinates.items() if v}
        if result.status != 'CLOSED':
            state.unknown.extend(result.obligations[:1])
        state.steps = result.expansions
        state.algorithm1 = result
        return result.status == 'CLOSED'

    def pattern_derivative_zero_certificate(self, rule, reaction, commutator):
        """Optional exact locality argument for a single source pattern."""
        return None
