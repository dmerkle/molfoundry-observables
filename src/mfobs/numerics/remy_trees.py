"""Exact propagation of the generated Rémy event model."""
from collections import defaultdict
from fractions import Fraction
from math import factorial
from typing import Iterable

from mfobs.models.remy.events import derive_jump_classes

def joint_distributions(
    generations: Iterable[int] = (3, 4, 5, 10, 50, 100),
) -> dict:
    """Integrate the derived weighted-growth operator on (P2,P3)."""

    requested = tuple(sorted(set(int(value) for value in generations)))
    if not requested or requested[0] < 0:
        raise ValueError("generations must be nonnegative")
    classes = derive_jump_classes()["classes"]
    weights: dict[tuple[int, int], int] = {(0, 0): 1}
    selected = {}
    moments = []
    for generation in range(requested[-1] + 1):
        total = sum(weights.values())
        expected_total = factorial(2 * generation) // factorial(generation)
        if total != expected_total:
            raise AssertionError("weighted Rémy mass is inconsistent")
        mean2 = Fraction(sum(p2 * value for (p2, _), value in weights.items()), total)
        mean3 = Fraction(sum(p3 * value for (_, p3), value in weights.items()), total)
        cross = Fraction(
            sum(p2 * p3 * value for (p2, p3), value in weights.items()), total
        )
        moments.append(
            {
                "generation": generation,
                "mean_P2": str(mean2),
                "mean_P3": str(mean3),
                "covariance_P2_P3": str(cross - mean2 * mean3),
                "mean_P2_float": float(mean2),
                "mean_P3_float": float(mean3),
            }
        )
        if generation in requested:
            selected[str(generation)] = {
                "generation": generation,
                "total_weight": str(total),
                "states": [
                    {
                        "P2": p2,
                        "P3": p3,
                        "weight": str(value),
                        "probability": value / total,
                    }
                    for (p2, p3), value in sorted(weights.items())
                ],
            }
        if generation == requested[-1]:
            break
        values_prefix = (2 * generation + 1, generation)
        next_weights: defaultdict[tuple[int, int], int] = defaultdict(int)
        for (p2, p3), weight in weights.items():
            values = (*values_prefix, p2, p3)
            for event in classes:
                propensity = sum(int(c) * v for c, v in zip(event["propensity"], values))
                if propensity < 0:
                    raise AssertionError("derived propensity became negative")
                if propensity:
                    delta = event["delta"]
                    next_weights[(p2 + delta[2], p3 + delta[3])] += (
                        weight * propensity
                    )
        weights = dict(next_weights)
    return {
        "status": "PASS",
        "method": "exact integer propagation of the derived jump-closure operator",
        "normalisation": "sum weights at generation n = (2n)!/n!",
        "generations": selected,
        "moments": moments,
    }
