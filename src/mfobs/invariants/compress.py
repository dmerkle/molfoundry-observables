"""Exact post-compression of an already closed, rate-separated linear system.

This is standard common-invariant-subspace/Krylov row-space computation, not
new graph rewriting or a general closure decision procedure. The graph-count
basis is replaced by rational linear combinations. Minimality as *functions*
is conditional on the source coordinates being independent functions on the
declared state language. No numerical parameter specialization is used.
"""
from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction
import hashlib
import json

import sympy as sp

from mfobs.core.linear_closure import RationalSpan as Span
from mfobs.core.exact_coefficients import rational_vector as clean


class CompressionError(ValueError):
    """The source is not an admissible exact closed linear-system result."""


def _rational(value, context):
    # Do not silently promote binary floating-point approximations to an
    # allegedly exact scientific certificate. Serialized decimal strings are
    # rational and are accepted; rate-dependent expressions are not.
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction, sp.Rational)):
        raise CompressionError(f'{context}: expected an exact rational, not {value!r}')
    try:
        return Fraction(str(value))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise CompressionError(f'{context}: expected an exact rational, got {value!r}') from exc


def _json_vector(row, names):
    return {names[i]: str(c) for i, c in sorted(row.items()) if c}


def _multiply(row, matrix):
    out = defaultdict(Fraction)
    for i, coefficient in row.items():
        for j, value in matrix[i].items():
            out[j] += coefficient * value
    return clean(out)


def _inputs(result):
    if not isinstance(result, dict) or result.get('closed') is not True:
        raise CompressionError('Compression requires a certified closed source result; UNKNOWN is not compressible')
    if result.get('unknown') or result.get('pending'):
        raise CompressionError('A closed source result must not have unknown or pending work')
    names = result.get('basis')
    if not isinstance(names, list) or not names or any(not isinstance(n, str) or not n for n in names) or len(set(names)) != len(names):
        raise CompressionError('basis must be a nonempty list of distinct observable names')
    index = {name: i for i, name in enumerate(names)}
    seeds = result.get('initial_observables')
    if seeds is None:
        depths = result.get('depths', {})
        seeds = [name for name in names if depths.get(name) == 0]
    if not isinstance(seeds, list) or not seeds or any(name not in index for name in seeds):
        raise CompressionError('Provide nonempty initial_observables, or source depths identifying the depth-zero seeds')
    seeds = list(dict.fromkeys(seeds))
    channels = result.get('rates')
    if not isinstance(channels, list):
        raise CompressionError('rates must be a list of channel declarations')
    by_id, matrices, rate_channels = {}, {}, {}
    for channel in channels:
        if not isinstance(channel, dict) or not isinstance(channel.get('id'), str) or not isinstance(channel.get('rate'), str) or not channel['rate']:
            raise CompressionError('Each rate channel requires string id and a nonempty rate-symbol name')
        rid, symbol = channel['id'], channel['rate']
        if rid in by_id:
            raise CompressionError(f'Duplicate rate channel id {rid!r}')
        by_id[rid] = symbol
        if symbol not in matrices:
            matrices[symbol] = [defaultdict(Fraction) for _ in names]
            rate_channels[symbol] = []
        rate_channels[symbol].append(rid)
    equations = result.get('equations_by_reaction')
    if not isinstance(equations, dict) or set(equations) != set(names):
        raise CompressionError('equations_by_reaction must explicitly include every source-basis row')
    for source, rows in equations.items():
        if not isinstance(rows, dict):
            raise CompressionError(f'Invalid equation row for {source}')
        for rid, values in rows.items():
            if rid not in by_id or not isinstance(values, dict):
                raise CompressionError(f'Unknown channel or malformed coefficients in {source}/{rid}')
            for target, coefficient in values.items():
                if target not in index:
                    raise CompressionError(f'Equation for {source} refers to unknown coordinate {target!r}')
                # Source equations already contain per_match multipliers, as
                # defined by outputs.matrix. Applying them again would be wrong.
                matrices[by_id[rid]][index[source]][index[target]] += _rational(coefficient, f'{source}/{rid}/{target}')
    matrices = {rate: [clean(row) for row in rows] for rate, rows in matrices.items()}
    return list(names), index, seeds, matrices, rate_channels


def compress(result):
    """Return a new exact linear-combination result; never mutate ``result``.

    For raw source coordinates f and rate-separated generator matrices A_r,
    discover the least rational row space W containing the supplied seed unit
    rows and closed under w -> w A_r for every independent rate symbol r.
    Channels with the same symbol are aggregated *before* space discovery.
    The returned transformation B has independent rows and satisfies
    B A_r = C_r B. Seed counts are recovered by the explicit seed_reconstruction.
    """
    names, index, seeds, matrices, rate_channels = _inputs(result)
    span = Span()
    rows, knames, origins, depths = [], [], {}, {}
    queue, steps = deque(), []

    def insert(candidate, origin, depth):
        residual, eliminated = span.reduce(candidate)
        if not residual:
            return None
        pivot = min(residual)
        scale = residual[pivot]
        row = {i: value / scale for i, value in residual.items()}
        name = f'K{len(rows)}'
        if not span.insert(name, row):
            raise AssertionError('Internal error: a nonzero exact residual must increase rank')
        rows.append(row)
        knames.append(name)
        queue.append(len(rows)-1)
        origins[name] = {
            **origin, 'candidate_in_source_coordinates': _json_vector(candidate, names),
            'eliminated_compressed_coordinates': {k: str(v) for k, v in eliminated.items()},
            'residual_scale': str(scale), 'pivot_source_coordinate': names[pivot],
        }
        depths[name] = depth
        return name

    for seed in seeds:
        insert({index[seed]: Fraction(1)}, {'kind': 'seed', 'seed': seed}, 0)
    while queue:
        i = queue.popleft()
        before = list(knames)
        for rate, matrix in matrices.items():
            insert(_multiply(rows[i], matrix),
                   {'kind': 'rate-generator-action', 'source': knames[i], 'rate': rate},
                   depths[knames[i]] + 1)
        steps.append({'expanded': knames[i], 'before': before, 'after': list(knames),
                      'queue': [knames[j] for j in queue]})
    equations = {}
    for name, row in zip(knames, rows):
        equations[name] = {}
        for rate, matrix in matrices.items():
            residual, coefficients = span.reduce(_multiply(row, matrix))
            if residual:
                raise AssertionError('Internal error: returned row space is not invariant')
            equations[name][rate] = {k: str(v) for k, v in coefficients.items() if v}
    seed_reconstruction = {}
    for seed in seeds:
        residual, coefficients = span.reduce({index[seed]: Fraction(1)})
        if residual:
            raise AssertionError('Internal error: a requested seed is missing')
        seed_reconstruction[seed] = {k: str(v) for k, v in coefficients.items() if v}
    initial_values = {}
    source_initial = result.get('initial_values') if result.get('initial_values_available', True) else None
    if source_initial is not None:
        if not isinstance(source_initial, dict) or set(source_initial) != set(names):
            raise CompressionError('When provided, source initial_values must cover exactly the source basis')
        initial = [_rational(source_initial[name], f'initial_values/{name}') for name in names]
        initial_values = {name: str(sum((value * initial[i] for i, value in row.items()), Fraction(0)))
                          for name, row in zip(knames, rows)}
    combinations = {name: _json_vector(row, names) for name, row in zip(knames, rows)}
    source_snapshot = json.dumps(result, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    out = {
        'schema_version': 1, 'model': str(result.get('model', 'linear-system')) + '-linear-compression',
        'title': str(result.get('title', result.get('model', 'Linear system'))) + ': exact rate-family linear compression',
        'closed': True, 'status': 'EXACT_RATE_FAMILY_LINEAR_COMPRESSION',
        'coordinate_kind': 'linear-combination-of-source-observables',
        'basis': knames, 'basis_combinations': combinations, 'basis_origins': origins,
        'initial_observables': [name for name in knames if depths[name] == 0],
        'source_initial_observables': seeds, 'seed_reconstruction': seed_reconstruction,
        'source_basis': names, 'source_dimension': len(names), 'dimension': len(knames),
        'dimension_reduction': len(names)-len(knames), 'depths': depths, 'steps': steps,
        'pending': [], 'unknown': [], 'initial_values': initial_values,
        'initial_values_available': source_initial is not None,
        'rates': [{'id': rate, 'rate': rate, 'per_match': '1', 'kind': 'abstract'} for rate in matrices],
        'aggregated_rate_channels': rate_channels, 'equations_by_reaction': equations,
        'normalisation': {name: {'scale': '1', 'meaning': ' + '.join(f'({coefficient})*{source}' for source, coefficient in combo.items())}
                          for name, combo in combinations.items()},
        'source_result_sha256': hashlib.sha256(source_snapshot.encode()).hexdigest(),
        'source_status': result.get('status'), 'source_proof_status': result.get('proof_status'),
        'proof_status': 'Exact rational row-worklist closure under every independent rate-symbol matrix after same-rate channel aggregation. The transformation and reconstructed seeds are checked exactly. Minimality as a function space is conditional on the source basis being functionally independent on its declared invariant state language; source closure and invariant preservation are inherited assumptions.',
        'minimality_scope': 'Smallest common invariant rational FUNCTION subspace containing the requested seeds, assuming independent original coordinate functions and treating distinct rate-symbol names as independently variable. Not minimal among individual graph-pattern bases, not a parameter-specialized rank claim, and not a general graph-closure decision.',
        'certificate': {
            'arithmetic': 'rational', 'rank': len(rows), 'source_coordinate_independence': 'assumed from source result',
            'invariance_identity': 'B A_rate = C_rate B for every aggregated rate symbol',
            'seed_identity': 'each source seed row = seed_reconstruction row times B',
            'construction': 'basis_origins records seed or rate action, prior-coordinate subtraction and nonzero residual scale; hence every added row belongs to the generated invariant space',
            'invariance_checked': True, 'seed_reconstruction_checked': True,
        },
        'raw_compositions': 0, 'ordered_composition_calls': 0, 'expansions': [],
        'limitations': [
            'post-processing requires an already closed finite exact source system',
            'compressed coordinates are linear combinations, not single graph patterns',
            'function-space minimality assumes the source coordinates are independent on the declared state language',
            'distinct rate-symbol names are treated independently; additional rate relations or numerical specializations may allow a smaller space',
            'no automatic proof of the source invariant language or source composition correctness is added',
        ],
    }
    # Independent dense checks make the sparse worklist certificate executable.
    transform = sp.Matrix([[sp.Rational(str(row.get(j, Fraction(0))))
                            for j in range(len(names))] for row in rows])
    if transform.rank() != len(rows):
        raise AssertionError('Internal error: dependent compressed coordinates')
    for rate, mat in matrices.items():
        original = sp.Matrix([[sp.Rational(str(row.get(j, Fraction(0))))
                               for j in range(len(names))] for row in mat])
        compressed = sp.Matrix([[sp.Rational(equations[name][rate].get(target, '0'))
                                 for target in knames] for name in knames])
        if transform * original != compressed * transform:
            raise AssertionError('Internal error: dense invariance check failed')
    for seed, coefficients in seed_reconstruction.items():
        reconstructed = sp.Matrix([[sp.Rational(coefficients.get(name, '0')) for name in knames]])
        expected = sp.zeros(1, len(names))
        expected[0, index[seed]] = 1
        if reconstructed * transform != expected:
            raise AssertionError('Internal error: dense seed reconstruction check failed')
    out['transformation_matrix'] = [[str(v) for v in transform.row(i)] for i in range(transform.rows)]
    from mfobs.outputs import matrix
    out['matrix'] = [[str(v) for v in row] for row in matrix(out).tolist()]
    return out
