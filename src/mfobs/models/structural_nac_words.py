"""Exact model operations."""
from __future__ import annotations
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product

from mfobs.core.semantic_contract import SemanticContract, theorem_preservation
from mfobs.inputs import InputError
from mfobs.models import structural_nac_oscillator as chemistry


def word_name(word: tuple[int, ...]) -> str:
    return "Q" + "".join(str(value) for value in word)


def word_meaning(word: tuple[int, ...]) -> str:
    labels = " → ".join(f"Q{value}" for value in word)
    return f"number of oriented consecutive occurrences of {labels}"


def concrete_rewrites(
    reaction: str,
) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """Expand a conditional graph rule into concrete word rewrites."""
    if reaction in chemistry.SWAP_TRANSITIONS:
        before, after = chemistry.SWAP_TRANSITIONS[reaction]
        return ((before, after),)
    attacker, old, new, guard = chemistry.RULE_TRANSITIONS[reaction]
    if guard is None:
        return (((attacker, old), (attacker, new)),)
    return tuple(
        ((attacker, old, probe), (attacker, new, probe))
        for probe in range(3) if probe != guard
    )


@lru_cache(maxsize=None)
def growth_generator_row(
    word: tuple[int, ...], reaction: str
) -> dict[tuple[int, ...], int]:
    """Exact count generator for inserting a copy immediately after each Qi site."""
    word = tuple(word)
    copied = chemistry.GROW_STATES[reaction]
    if len(word) == 1:
        return {(copied,): 1} if word == (copied,) else {}

    # An insertion after the central anchor replaces L-1 old windows by L new
    # windows. Their difference depends only on this 2L-2 site context. A
    # nonzero row can occur only when at least one of those old/new windows is
    # the queried word, so construct that union of matching cylinders directly
    # rather than enumerating all 3^(2L-3) contexts (almost all of which cancel).
    length = len(word)
    context_length = 2 * length - 2
    anchor = length - 2
    candidates: set[tuple[int, ...]] = set()

    def add_matching_contexts(start: int, *, after_insertion: bool) -> None:
        fixed = {anchor: copied}
        for offset, value in enumerate(word):
            index = start + offset
            if after_insertion and index == anchor + 1:
                if value != copied:
                    return
                continue
            before_index = index - 1 if after_insertion and index > anchor + 1 else index
            if before_index in fixed and fixed[before_index] != value:
                return
            fixed[before_index] = value
        free = [index for index in range(context_length) if index not in fixed]
        for values in product(range(3), repeat=len(free)):
            context = [0] * context_length
            for index, value in fixed.items():
                context[index] = value
            for index, value in zip(free, values):
                context[index] = value
            candidates.add(tuple(context))

    for start in range(length - 1):
        add_matching_contexts(start, after_insertion=False)
    for start in range(length):
        add_matching_contexts(start, after_insertion=True)

    result: Counter[tuple[int, ...]] = Counter()
    for before in candidates:
        after = before[:anchor + 1] + (copied,) + before[anchor + 1:]
        loss = sum(
            before[start:start + length] == word for start in range(length - 1)
        )
        gain = sum(
            after[start:start + length] == word for start in range(length)
        )
        if gain != loss:
            result[before] += gain - loss
    return {
        target: coefficient
        for target, coefficient in sorted(result.items())
        if coefficient
    }


@lru_cache(maxsize=None)
def generator_row(word: tuple[int, ...], reaction: str) -> dict[tuple[int, ...], int]:
    """Exact generator image of one cyclic-word count in one rate channel.

    The changed site can occupy each position of ``word``.  For every such
    overlap, the union of the rule window and the observable window is
    enumerated symbolically.  Each union assignment is itself one connected
    graph observable, so the returned integer row is exact for every declared
    macrocycle longer than the union window.
    """
    word = tuple(word)
    if not word or any(value not in (0, 1, 2) for value in word):
        raise ValueError("word must be a nonempty Q0/Q1/Q2 tuple")
    if reaction not in chemistry.REACTION_ORDER:
        raise ValueError(f"unknown reaction {reaction!r}")
    if reaction in chemistry.GROW_STATES:
        return growth_generator_row(word, reaction)
    result: Counter[tuple[int, ...]] = Counter()
    for precondition, postcondition in concrete_rewrites(reaction):
        changed = [
            index for index, (before, after) in
            enumerate(zip(precondition, postcondition)) if before != after
        ]
        observable_starts = sorted({
            changed_site - observable_index
            for changed_site in changed
            for observable_index in range(len(word))
        })
        for observable_start in observable_starts:
            observable_stop = observable_start + len(word) - 1
            union_start = min(0, observable_start)
            union_stop = max(len(precondition) - 1, observable_stop)
            constrained = {index: value for index, value in enumerate(precondition)}
            free = [
                index for index in range(union_start, union_stop + 1)
                if index not in constrained
            ]
            for values in product(range(3), repeat=len(free)):
                before = dict(constrained)
                before.update(zip(free, values))
                before_window = tuple(
                    before[index]
                    for index in range(observable_start, observable_stop + 1)
                )
                after = dict(before)
                after.update(enumerate(postcondition))
                after_window = tuple(
                    after[index]
                    for index in range(observable_start, observable_stop + 1)
                )
                coefficient = int(after_window == word) - int(before_window == word)
                if coefficient:
                    union_word = tuple(
                        before[index] for index in range(union_start, union_stop + 1)
                    )
                    result[union_word] += coefficient
    return {
        target: coefficient
        for target, coefficient in sorted(result.items())
        if coefficient
    }


class Adapter:
    """Worklist adapter exposing exact connected graph-word observables."""

    semantic_contract_required = True

    def __init__(self, inputs):
        inputs.validate_overlap_policy()
        configuration = inputs.read("reduction")
        if set(configuration) != {"method", "coordinate_family"}:
            raise InputError("Unexpected structural-word reduction fields")
        if configuration["method"] != "exact-oriented-cyclic-word-generator-v1":
            raise InputError("Wrong structural-word reduction method")
        if configuration["coordinate_family"] != "all-connected-Q0-Q1-Q2-words":
            raise InputError("Wrong structural-word coordinate family")

        language = inputs.read("state_language")
        if language.get("kind") != "unbounded-oriented-covalent-switch-macrocycles":
            raise InputError("Wrong structural oscillator state language")
        if language.get("alphabet") != ["Q0", "Q1", "Q2"]:
            raise InputError("Structural oscillator alphabet must be Q0/Q1/Q2")
        minimum = language.get("minimum_ring_size")
        if type(minimum) is not int or minimum < 5:
            raise InputError("minimum_ring_size must be an integer at least five")

        expected = chemistry.rules()
        actual = inputs.engine_rules()
        if set(actual) != set(expected):
            raise InputError("Input must contain exactly the certified rule families")
        for name, rule in expected.items():
            expected_engine = rule
            if actual[name].canonKey != expected_engine.canonKey:
                raise InputError(f"Rule {name} differs from the certified DPO rule")
        supplied_factors = {
            item["id"]: Fraction(item["per_match"])
            for item in inputs.data["rules"]
        }
        if supplied_factors != chemistry.rate_factors():
            raise InputError("Per-match factors do not normalize hydrogen choices")

        self.inputs = inputs
        self.language = language
        self.words: dict[str, tuple[int, ...]] = {}
        self.names: dict[tuple[int, ...], str] = {}
        self.origins: dict[str, dict] = {}
        self.discoveries: list[dict] = []
        self.generator_rows: dict[str, dict] = {}

        scope = (
            f"all oriented macrocycles of size N >= {minimum}, whenever the "
            "counted word and its exact generator windows have length <= N"
        )
        preservation = []
        for name, (_attacker, old, new, guard) in chemistry.RULE_TRANSITIONS.items():
            preservation.append(
                theorem_preservation(
                    name,
                    scope,
                    {
                        "literal_rule": "canonical key equals the generated conditional DPO rule",
                        "scaffold": "all atoms and every C-S-O-C backbone edge are preserved",
                        "state_change": (
                            f"only the target chord changes Q{old}->Q{new}; exactly one "
                            "allowed chord remains at every site"
                        ),
                        "guard": (
                            "none" if guard is None else
                            f"NoEdge(e{guard}) on an unchanged probe is equivalent to probe != Q{guard}"
                        ),
                    },
                )
            )
        for name, (before, after) in chemistry.SWAP_TRANSITIONS.items():
            preservation.append(
                theorem_preservation(
                    name,
                    scope,
                    {
                        "literal_rule": "canonical key equals the generated conditional DPO rule",
                        "scaffold": "all atoms and every C-S-O-C backbone edge are preserved",
                        "state_change": f"adjacent chord states {before} are exchanged to {after}",
                        "guard": "none",
                    },
                )
            )
        for name, state in chemistry.GROW_STATES.items():
            preservation.append(
                theorem_preservation(
                    name,
                    scope,
                    {
                        "literal_rule": "canonical key equals the generated DPO condensation rule",
                        "atoms": "the capped precursor is consumed and its cap leaves as mapped water",
                        "state_change": f"one new Q{state} repeat is inserted after a Q{state} template site",
                        "size_change": "the macrocycle grows from N to N+1; no graph atom is created",
                        "reservoir": "precursor activity is chemostatted into the effective rate beta",
                    },
                )
            )
        self.semantic_contract = SemanticContract(
            domain=language["kind"],
            reducer="exact symbolic local-word generator on oriented cycles",
            exact_reduction_scope=(
                "Connected word observables whose rule/observable union fits in "
                "the macrocycle; integer coefficients, no moment factorization."
            ),
            preservation=tuple(preservation),
        )

    def _identity_rule(self, word: tuple[int, ...]):
        conditional = chemistry.observable_rule(word)
        return conditional

    def _add_word(self, word, state, *, origin, supplied_rule=None):
        word = tuple(word)
        if word in self.names:
            return self.names[word]
        name = word_name(word)
        if name in state.patterns:
            raise InputError(f"Structural word name collision: {name}")
        rule = supplied_rule if supplied_rule is not None else self._identity_rule(word)
        state.add(name, rule, word_meaning(word))
        self.words[name] = word
        self.names[word] = name
        self.origins[name] = origin
        return name

    def initialize(self, state, composer):
        observables = self.inputs.engine_observables()
        for name in self.inputs.data["initial_observables"]:
            if len(name) < 2 or name[0] != "Q" or any(value not in "012" for value in name[1:]):
                raise InputError(
                    "Cyclic-word observables use canonical names Q followed by "
                    "their Q0/Q1/Q2 word, for example Q0 or Q012"
                )
            word = tuple(int(value) for value in name[1:])
            expected = self._identity_rule(word)
            if observables[name].canonKey != expected.canonKey:
                raise InputError(
                    f"Seed {name} is not the literal graph identity for word {name[1:]}"
                )
            self._add_word(
                word,
                state,
                supplied_rule=observables[name],
                origin={"kind": "input-seed", "word": list(word)},
            )

    def run_closure(self, state, composer, *, max_expansions, maximum_dimension):
        from mfobs.core.linear_closure import close_linear_span
        from mfobs.core.exact_coefficients import add_vectors as add

        forms = {name: {word: Fraction(1)} for name, word in self.words.items()}
        channels = {}
        for reaction in self.inputs.data['rules']:
            channels.setdefault(reaction['rate'], []).append(reaction)

        def image(source, channel):
            state.current_depth = state.patterns[source].depth
            total = {}
            for reaction in channels[channel]:
                # generator_row counts normalized word events: the theorem
                # already removes mapped-hydrogen match multiplicities.
                for word, coefficient in forms[source].items():
                    total = add(total, generator_row(word, reaction['id']),
                                coefficient)
            return total, {'method': 'exact local-word generator, summed by rate channel'}

        def added(name, vector, source, channel, before, scale):
            forms[name] = vector
            singleton = next(iter(vector.items())) if len(vector) == 1 else None
            rule = self._identity_rule(singleton[0]) if singleton and singleton[1] == 1 else None
            state.add(name, rule, 'exact linear combination of cyclic word counts')

        result = close_linear_span(
            seeds=forms.copy(), channels=channels, image=image,
            encode=lambda v: {word_name(w): str(c) for w, c in v.items()},
            designate=lambda v, i: f'O{i}', on_add=added,
            on_expansion=state.record_expansion,
            maximum_dimension=maximum_dimension, maximum_expansions=max_expansions,
        )
        self.word_forms = forms
        self.origins = result.origins
        self.discoveries = result.steps
        self.generator_rows = result.steps
        state.algorithm1 = result
        state.obligations = result.obligations
        state.queue.clear()
        state.queue.extend(result.pending)
        state.channel_equations = {}
        for step in result.steps:
            if 'coordinates' in step:
                state.channel_equations.setdefault(step['source'], {})[step['channel']] = step['coordinates']
        state.steps = result.expansions
        if result.status != 'CLOSED':
            state.unknown.extend(result.obligations[:1])
        # Retain literal operands for inspection, without inventing one graph
        # representative for a linear combination.
        rules = {w: self._identity_rule(w) for v in forms.values() for w in v}
        self.expressions = {n: {word_name(w): (c, rules[w]) for w, c in v.items()}
                            for n, v in forms.items()}
        return result.status == 'CLOSED'

    def finalize(self, state, composer, complete):
        initial = self.inputs.data["initial_state"]
        if set(initial) != {"word"}:
            raise InputError("initial_state must contain exactly one cyclic word")
        initial_word = initial["word"]
        if (
            not isinstance(initial_word, list)
            or len(initial_word) < self.language["minimum_ring_size"]
            or any(type(value) is not int or value not in (0, 1, 2)
                   for value in initial_word)
        ):
            raise InputError("initial cyclic word is outside the declared language")
        initial_word = tuple(initial_word)
        coordinates, values = {}, {}
        for name, form in self.word_forms.items():
            coordinates[name] = {'terms': {word_name(w): str(c) for w, c in form.items()},
                                 'meaning': 'exact linear combination of cyclic word counts'}
            values[name] = str(sum(c * chemistry.word_count(initial_word, w) for w, c in form.items()))
        maximum = max((len(w) for v in self.word_forms.values() for w in v), default=0)
        return {
            "status": "CLOSED_EXACT_WORD_SYSTEM" if complete else "UNKNOWN",
            "algorithm1_status": state.algorithm1.status,
            "proof_status": (
                "Every row is the exact conditional-rule generator on the "
                "declared cyclic-word language. No factorisation is used."
            ),
            "word_coordinates": coordinates,
            "basis_origins": self.origins,
            "initial_values": values,
            "maximum_generated_word_length": maximum,
            "synthesis_state_evaluations": 0,
            "synthesis_graph_compositions": 0,
            "composition_note": (
                "The literal GML rules are canonical-key checked. Their exact "
                "action is then reduced by the supplied infinite cyclic-word "
                "theorem, rather than enumerating finite host graphs."
            ),
            "limitations": [
                "The material, net ring shifts, and precursor condensation are a chemically typed hypothetical benchmark.",
                "Unbounded growth assumes chemostatted activity of the three explicit capped precursors.",
                "A budget-limited run is UNKNOWN; finite hierarchy integrations are explicit approximations.",
                "The exact word theorem applies only when each local union window fits inside N.",
            ],
        }
