"""Reusable incremental exact linear span over rational vectors."""
from __future__ import annotations

from fractions import Fraction
from typing import Sequence


Vector = tuple[Fraction, ...]


class ExactIncrementalSpan:
    """Maintain an exact RREF basis and insertion-order coordinates."""

    def __init__(self, first: Sequence[int | Fraction]) -> None:
        vector = tuple(Fraction(value) for value in first)
        if not any(vector):
            raise ValueError("the first span vector must be nonzero")
        pivot = next(index for index, value in enumerate(vector) if value)
        scale = vector[pivot]
        self.vector_dimension = len(vector)
        self.pivots = [pivot]
        self.rows = [[value / scale for value in vector]]
        self.representations = [[Fraction(1, 1) / scale]]

    @property
    def dimension(self) -> int:
        return len(self.rows)

    def _reduce(
        self, target: Sequence[int | Fraction]
    ) -> tuple[list[Fraction], list[Fraction]]:
        if len(target) != self.vector_dimension:
            raise ValueError("span vector has the wrong dimension")
        residual = [Fraction(value) for value in target]
        coordinates = [Fraction() for _ in self.rows]
        for pivot, row, representation in zip(
            self.pivots, self.rows, self.representations
        ):
            factor = residual[pivot]
            if not factor:
                continue
            residual = [value - factor * basis_value
                        for value, basis_value in zip(residual, row)]
            coordinates = [value + factor * basis_value
                           for value, basis_value in zip(coordinates, representation)]
        return residual, coordinates

    def solve(
        self, target: Sequence[int | Fraction]
    ) -> tuple[Fraction, ...] | None:
        residual, coordinates = self._reduce(target)
        return tuple(coordinates) if not any(residual) else None

    def add(self, target: Sequence[int | Fraction]) -> None:
        """Add one vector already known to be independent."""

        residual, old_coordinates = self._reduce(target)
        if not any(residual):
            raise ValueError("cannot add a dependent span vector")
        new_pivot = next(index for index, value in enumerate(residual) if value)
        scale = residual[new_pivot]
        new_row = [value / scale for value in residual]
        for representation in self.representations:
            representation.append(Fraction())
        new_representation = [
            -value / scale for value in old_coordinates
        ] + [Fraction(1, 1) / scale]

        for index, row in enumerate(self.rows):
            factor = row[new_pivot]
            if factor:
                self.rows[index] = [
                    value - factor * new_value
                    for value, new_value in zip(row, new_row)
                ]
                self.representations[index] = [
                    value - factor * new_value
                    for value, new_value in zip(
                        self.representations[index], new_representation
                    )
                ]

        insertion = next(
            (index for index, pivot in enumerate(self.pivots) if pivot > new_pivot),
            len(self.pivots),
        )
        self.pivots.insert(insertion, new_pivot)
        self.rows.insert(insertion, new_row)
        self.representations.insert(insertion, new_representation)
