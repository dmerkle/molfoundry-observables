"""Algorithm 1: exact, channel-wise residual closure.

The image oracle sums *all* rules of a rate channel, including their per-match
multipliers. The span and normalisation implement the backend's exact algebra;
neither may use sampling to certify membership. One nonzero residual adds one
observable, even when that observable is a linear combination of patterns.
"""
from mfobs.core.exact_coefficients import add_vectors

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd, lcm


class UnsupportedImage(Exception):
    """The backend cannot certify this generator image."""


class RationalSpan:
    """Sparse exact elimination with coordinates in the inserted basis."""

    def __init__(self):
        self.rows = []

    add = staticmethod(add_vectors)

    def reduce(self, vector):
        residual = {key: Fraction(value) for key, value in vector.items() if value}
        coordinates = {}
        for pivot, row, expression in self.rows:
            value = residual.get(pivot, 0)
            if value:
                residual = self.add(residual, row, -value)
                coordinates = self.add(coordinates, expression, value)
        return residual, coordinates

    def insert(self, name, vector):
        residual, coordinates = self.reduce(vector)
        if not residual:
            return False
        pivot = min(residual)
        value = residual[pivot]
        self.rows.append((pivot,
                          {key: c / value for key, c in residual.items()},
                          {**{key: -c / value for key, c in coordinates.items()},
                           name: Fraction(1) / value}))
        return True


def primitive(vector):
    """Return a primitive rational vector and its exact original scale."""
    denominator = lcm(*(Fraction(c).denominator for c in vector.values()))
    divisor = reduce(gcd, (abs(int(c * denominator)) for c in vector.values()))
    scale = Fraction(divisor, denominator)
    if vector[min(vector)] < 0:
        scale = -scale
    return {key: value / scale for key, value in vector.items()}, scale


@dataclass
class LinearClosure:
    status: str
    basis: dict
    origins: dict
    steps: list
    pending: list
    obligations: list
    expansions: list


def close_linear_span(*, seeds, channels, image, span=None, normalize=primitive,
                      encode=lambda v: {str(k): str(c) for k, c in v.items()},
                      designate=lambda v, i: f'O{i}', on_add=None, on_step=None,
                      on_expansion=None,
                      maximum_dimension=64, maximum_expansions=None):
    """Execute the paper's FIFO loop; callbacks provide algebra and evidence.

    ``image(name, channel)`` returns a sparse exact vector and trace metadata.
    ``on_add`` attaches a backend representation to the residual, before it is
    expanded. ``normalize`` may choose a simpler representative modulo the
    current span; the core verifies its residual is exactly the same direction.
    Seeds are kept verbatim, including scale. Budget exits retain the
    failed image plus every unprocessed (source, channel) obligation.
    An empty channel set represents the zero generator.
    """
    channels = tuple(channels)
    if not seeds or any(not c for c in channels):
        raise ValueError('nonempty seeds and nonempty rate-channel names are required')
    if maximum_dimension < len(seeds):
        raise ValueError('maximum_dimension is smaller than the seed basis')
    if maximum_expansions is not None and maximum_expansions < 0:
        raise ValueError('maximum_expansions must be nonnegative')
    span = span if span is not None else RationalSpan()
    basis, origins = {}, {}
    for name, vector in seeds.items():
        if not span.insert(name, vector):
            raise ValueError('seed observables must be nonzero and linearly independent')
        basis[name] = vector
        origins[name] = {'origin': 'input-seed'}
    queue = deque(basis)
    steps = []
    expanded = 0
    expansions = []

    def stopped(status, source, index, detail):
        obligations = [{'source': source, 'channel': c} for c in channels[index:]]
        obligations += [{'source': n, 'channel': c} for n in queue for c in channels]
        obligations[0].update(detail)
        return LinearClosure(status, basis, origins, steps,
                             [source, *queue], obligations, expansions)

    while queue:
        source = queue.popleft()
        basis_before = list(basis)
        if channels and maximum_expansions is not None and expanded >= maximum_expansions:
            return stopped('UNKNOWN_EXPANSION_LIMIT', source, 0,
                           {'reason': 'expansion budget exhausted'})
        for index, channel in enumerate(channels):
            try:
                vector, evidence = image(source, channel)
            except UnsupportedImage as exc:
                return stopped('UNKNOWN_REDUCTION', source, index, {'reason': str(exc)})
            residual, before = span.reduce(vector)
            step = {'source': source, 'channel': channel, **evidence,
                    'image': encode(vector),
                    'coordinates_before_extension': {n: str(c) for n, c in before.items()},
                    'residual_before_extension': encode(residual),
                    'added': None, 'residual_scale': None}
            if residual:
                if len(basis) >= maximum_dimension:
                    step['status'] = 'UNKNOWN_DIMENSION_LIMIT'
                    steps.append(step)
                    return stopped('UNKNOWN_DIMENSION_LIMIT', source, index,
                                   {'reason': 'dimension budget exhausted',
                                    'image': encode(vector), 'residual': encode(residual)})
                discovered, scale = normalize(residual)
                # A backend may choose a simpler representative modulo the
                # current span, but MUST add precisely this one new direction.
                representative_residual, representative_coordinates = span.reduce(discovered)
                difference, _ = span.reduce({
                    k: residual.get(k, 0) - scale * representative_residual.get(k, 0)
                    for k in residual.keys() | representative_residual.keys()})
                if any(difference.values()):
                    raise AssertionError('extension does not represent exactly the residual direction')
                step['representative_coordinates'] = {n: str(c) for n, c in representative_coordinates.items()}
                name = designate(discovered, len(basis))
                if name in basis:
                    raise ValueError(f'duplicate discovered observable name: {name}')
                if not span.insert(name, discovered):
                    raise AssertionError('residual failed to extend the exact span')
                basis[name] = discovered
                origins[name] = {'origin': 'computed-channel-residual',
                                 'discovered_from': source, 'rate_channel': channel}
                if on_add is not None:
                    on_add(name, discovered, source, channel, before, scale)
                queue.append(name)
                step.update(added=name, residual_scale=str(scale))
            remainder, coordinates = span.reduce(vector)
            if remainder:
                raise AssertionError('image remains outside the extended span')
            step['coordinates'] = {n: str(c) for n, c in coordinates.items() if c}
            steps.append(step)
            if on_step is not None:
                on_step(step)
        expansions.append({'expanded': source, 'before': basis_before,
                           'after': list(basis), 'queue': list(queue)})
        if on_expansion is not None:
            on_expansion(expansions[-1])
        expanded += 1
    return LinearClosure('CLOSED', basis, origins, steps, [], [], expansions)
