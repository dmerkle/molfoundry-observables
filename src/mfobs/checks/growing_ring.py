"""Concrete checks of the growing-ring word generator."""
from __future__ import annotations
from mfobs.models import structural_nac_oscillator as chemistry


def evaluate_row(
    row: dict[tuple[int, ...], int], states: tuple[int, ...]
) -> int:
    return sum(
        coefficient * chemistry.word_count(states, target)
        for target, coefficient in row.items()
    )


def direct_jump_sum(
    word: tuple[int, ...], reaction: str, states: tuple[int, ...]
) -> int:
    before = chemistry.word_count(states, word)
    return sum(
        chemistry.word_count(target, word) - before
        for target in chemistry.transition_events(states, reaction)
    )
