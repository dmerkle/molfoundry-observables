"""Generic terminal-growth word grammar for sequence observables.

This is the reusable combinatorial bridge between recursive graph grammars and
polymer microstructure.  A caller supplies an alphabet and motif words.  The
module constructs the smallest fixed suffix depth needed to classify every
motif jump, augments the observable coordinates with exact suffix indicators,
and emits a certified finite event partition for words of arbitrary length.

No chemical names, equations, closure dimensions, or finite word lists are
inputs to the construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Mapping, Sequence

from mfobs.core.polynomial_jump_closure import (
    JumpEvent,
    JumpModelProvenance,
    PolynomialJumpModel,
    PolynomialRing,
)


@dataclass(frozen=True)
class WordObservable:
    name: str
    kind: str
    word: tuple[str, ...] = ()
    symbol: str | None = None
    role: str = "selected"
    meaning: str = ""


@dataclass(frozen=True)
class WordState:
    word: tuple[str, ...]
    pools: tuple[int, ...]
    gate_open: bool = True


@dataclass(frozen=True)
class TerminalGrowthGrammar:
    """Open-system terminal polymer growth with an optional Boolean gate."""

    name: str
    alphabet: tuple[str, ...]
    motifs: tuple[tuple[str, ...], ...]
    gated: bool
    input_definition_sha256: str

    def __post_init__(self) -> None:
        if not self.alphabet or len(self.alphabet) != len(set(self.alphabet)):
            raise ValueError("alphabet symbols must be nonempty and unique")
        alphabet = set(self.alphabet)
        if any(len(motif) < 2 or set(motif) - alphabet for motif in self.motifs):
            raise ValueError("motifs must be alphabet words of length at least two")
        if len(self.motifs) != len(set(self.motifs)):
            raise ValueError("motifs must be unique")

    @property
    def suffix_depth(self) -> int:
        return max((len(motif) - 1 for motif in self.motifs), default=0)

    @property
    def suffix_classes(self) -> tuple[tuple[str, ...], ...]:
        depth = self.suffix_depth
        if not depth:
            return ((),)
        short = [
            word
            for length in range(depth)
            for word in product(self.alphabet, repeat=length)
        ]
        full = list(product(self.alphabet, repeat=depth))
        return tuple(short + full)

    @staticmethod
    def suffix_name(word: tuple[str, ...]) -> str:
        return "S_" + ("empty" if not word else "".join(word))

    def classify(self, word: Sequence[str]) -> tuple[str, ...]:
        word = tuple(word)
        if not self.suffix_depth:
            return ()
        if len(word) < self.suffix_depth:
            return word
        return word[-self.suffix_depth :]

    def primary_observables(self) -> tuple[WordObservable, ...]:
        result = [
            WordObservable(
                f"N_{symbol}", "symbol-count", symbol=symbol,
                meaning=f"number of incorporated {symbol} repeat units",
            )
            for symbol in self.alphabet
        ]
        result.extend(
            WordObservable(
                "C_" + "".join(motif), "factor-count", word=motif,
                meaning=f"number of consecutive {''.join(motif)} sequence motifs",
            )
            for motif in self.motifs
        )
        result.extend(
            WordObservable(
                f"M_{symbol}", "pool-count", symbol=symbol,
                role="support",
                meaning=f"number of free {symbol} monomers",
            )
            for symbol in self.alphabet
        )
        if self.gated:
            result.append(WordObservable(
                "X", "gate-open", role="support",
                meaning="indicator that the hemilabile donor is unbound"
            ))
        return tuple(result)

    def auxiliary_observables(self) -> tuple[WordObservable, ...]:
        if not self.suffix_depth:
            return ()
        # The empty class is represented by 1-sum(nonempty classes), keeping
        # the indicator coordinate family linearly independent.
        return tuple(
            WordObservable(
                self.suffix_name(word), "suffix-class", word=word,
                role="auxiliary",
                meaning=(
                    f"indicator that the current word ends in {''.join(word)}; "
                    f"generated from motif depth {self.suffix_depth + 1}"
                ),
            )
            for word in self.suffix_classes if word
        )

    def observables(self) -> tuple[WordObservable, ...]:
        return self.primary_observables() + self.auxiliary_observables()

    def coordinate_values(self, state: WordState) -> dict[str, int]:
        if len(state.pools) != len(self.alphabet) or any(value < 0 for value in state.pools):
            raise ValueError("pool inventory differs from the grammar alphabet")
        if set(state.word) - set(self.alphabet):
            raise ValueError("state word contains a symbol outside the grammar alphabet")
        result: dict[str, int] = {}
        for observable in self.observables():
            if observable.kind == "symbol-count":
                result[observable.name] = state.word.count(observable.symbol)
            elif observable.kind == "factor-count":
                width = len(observable.word)
                result[observable.name] = sum(
                    state.word[index:index + width] == observable.word
                    for index in range(len(state.word) - width + 1)
                )
            elif observable.kind == "pool-count":
                result[observable.name] = state.pools[self.alphabet.index(observable.symbol)]
            elif observable.kind == "gate-open":
                result[observable.name] = int(state.gate_open)
            elif observable.kind == "suffix-class":
                width = len(observable.word)
                result[observable.name] = int(
                    len(state.word) >= width
                    and state.word[-width:] == observable.word
                )
            else:
                raise AssertionError(f"unknown word observable {observable.kind!r}")
        return result

    def suffix_indicator(self, ring: PolynomialRing, suffix: tuple[str, ...]):
        """Return the exact indicator of one mechanically generated suffix class."""

        if suffix not in self.suffix_classes:
            raise ValueError("suffix is not a class of this grammar")
        if not self.suffix_depth:
            return ring.constant(1)
        if len(suffix) == self.suffix_depth:
            return ring.coordinate(self.suffix_name(suffix))
        if suffix:
            value = ring.coordinate(self.suffix_name(suffix))
            extensions = [
                ring.coordinate(self.suffix_name((symbol,) + suffix))
                for symbol in self.alphabet
            ]
            return ring.add(value, *(ring.scale(item, -1) for item in extensions))
        length_one = [
            ring.coordinate(self.suffix_name((symbol,)))
            for symbol in self.alphabet
        ]
        return ring.add(ring.constant(1), *(ring.scale(value, -1) for value in length_one))

    def _jump(
        self,
        source_suffix: tuple[str, ...],
        *,
        appended: str | None = None,
        pool_delta: Mapping[str, int] | None = None,
        gate_delta: int = 0,
    ) -> tuple[int, ...]:
        old_word = source_suffix
        new_word = old_word + ((appended,) if appended is not None else ())
        old = WordState(old_word, (2,) * len(self.alphabet), self.gated)
        new_pools = tuple(
            value + (pool_delta or {}).get(symbol, 0)
            for symbol, value in zip(self.alphabet, old.pools)
        )
        new = WordState(new_word, new_pools, bool(int(old.gate_open) + gate_delta))
        before = self.coordinate_values(old)
        after = self.coordinate_values(new)
        return tuple(after[name] - before[name] for name in before)

    def polynomial_ring(self) -> PolynomialRing:
        """Construct the finite coordinate algebra, independently of events."""

        observables = self.observables()
        names = tuple(observable.name for observable in observables)
        suffix_names = tuple(
            observable.name for observable in observables
            if observable.kind == "suffix-class"
        )
        boolean = suffix_names + (("X",) if self.gated else ())
        groups = tuple(
            tuple(self.suffix_name(word) for word in self.suffix_classes if len(word) == length)
            for length in range(1, self.suffix_depth + 1)
        )
        suffix_words = tuple(word for word in self.suffix_classes if word)
        products = {}
        for index, left in enumerate(suffix_words):
            for right in suffix_words[index + 1:]:
                pair = (self.suffix_name(left), self.suffix_name(right))
                if len(left) >= len(right) and left[-len(right):] == right:
                    products[pair] = self.suffix_name(left)
                elif len(right) >= len(left) and right[-len(left):] == left:
                    products[pair] = self.suffix_name(right)
                else:
                    products[pair] = None
        return PolynomialRing(
            names,
            boolean_variables=boolean,
            exclusive_groups=groups,
            boolean_products=products,
        )

    def polynomial_jump_model(self) -> PolynomialJumpModel:
        ring = self.polynomial_ring()
        names = ring.variables
        events: list[JumpEvent] = []

        for symbol in self.alphabet:
            delta = [0] * len(names)
            delta[names.index(f"M_{symbol}")] = 1
            events.append(JumpEvent(
                f"feed-{symbol}", f"beta_{symbol}", tuple(delta),
                ring.constant(1), f"feed_{symbol}", "one persistent reservoir match",
            ))

        for suffix in self.suffix_classes:
            indicator = self.suffix_indicator(ring, suffix)
            for symbol in self.alphabet:
                propensity = ring.multiply(ring.coordinate(f"M_{symbol}"), indicator)
                if self.gated:
                    propensity = ring.multiply(propensity, ring.coordinate("X"))
                events.append(JumpEvent(
                    f"append-{symbol}-from-{''.join(suffix) or 'empty'}",
                    f"k_{symbol}",
                    self._jump(suffix, appended=symbol, pool_delta={symbol: -1}),
                    propensity,
                    f"insert_{symbol}",
                    (
                        f"suffix={''.join(suffix) or 'empty'}"
                        + (" and NoEdge(Pd,N)" if self.gated else "")
                    ),
                ))

        if self.gated:
            delta = [0] * len(names)
            delta[names.index("X")] = -1
            events.append(JumpEvent(
                "donor-bind", "a", tuple(delta), ring.coordinate("X"),
                "donor_bind", "X=1 and NoEdge(Pd,N)",
            ))
            delta[names.index("X")] = 1
            events.append(JumpEvent(
                "donor-unbind", "b", tuple(delta),
                ring.subtract(ring.constant(1), ring.coordinate("X")),
                "donor_unbind", "X=0",
            ))

        certificate = {
            "status": "VERIFIED_ALL_STATES",
            "method": "finite suffix congruence induced mechanically by the longest selected motif",
            "alphabet": list(self.alphabet),
            "motifs": ["".join(motif) for motif in self.motifs],
            "suffix_depth": self.suffix_depth,
            "suffix_classes": ["".join(value) or "empty" for value in self.suffix_classes],
            "argument": (
                "Appending one symbol changes a motif of length at most k only through "
                "the previous suffix of length at most k-1. Pool and gate jumps are "
                "constant. Hence every word of every length belongs to exactly one "
                "listed event class and has the recorded jump."
            ),
            "enumeration_used": False,
        }
        return PolynomialJumpModel(
            self.name, ring, tuple(events), certificate,
            self.input_definition_sha256,
            JumpModelProvenance(
                backend="recursive terminal-word semantic oracle",
                event_partition_origin=(
                    "finite suffix congruence generated from alphabet and motifs"
                ),
                jump_origin="exact before/after evaluation of word coordinates",
                activity_origin="exact suffix, pool, and Boolean word algebra",
                jump_projection="semantic word transition",
                composition_scope="none; independent validation oracle",
                all_jumps_composition_derived=False,
                supplied_delta_table=False,
            ),
        )

    def successors(self, state: WordState) -> tuple[tuple[str, WordState], ...]:
        """Concrete semantic transitions, used only for report expansion/tests."""

        result = []
        for index, symbol in enumerate(self.alphabet):
            pools = list(state.pools)
            pools[index] += 1
            result.append((f"feed_{symbol}", WordState(state.word, tuple(pools), state.gate_open)))
            if state.pools[index] and (state.gate_open or not self.gated):
                pools = list(state.pools)
                pools[index] -= 1
                result.append((
                    f"insert_{symbol}",
                    WordState(state.word + (symbol,), tuple(pools), state.gate_open),
                ))
        if self.gated:
            result.append((
                "donor_bind" if state.gate_open else "donor_unbind",
                WordState(state.word, state.pools, not state.gate_open),
            ))
        return tuple(result)
