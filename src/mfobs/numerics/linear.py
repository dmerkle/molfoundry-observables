"""Integrate closed linear discovery records; no event partition is required."""
from fractions import Fraction
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from mfobs.outputs import matrix


def linear_matrix(record, rates):
    """Read exported coefficients, or channel-separated Algorithm 1 steps.

    For worklist steps, rate keys are channel names (e.g. G for Rémy's
    sum of six insertion rules). For exports, keys are parameter names.
    """
    if not record.get('closed', record.get('status') == 'CLOSED') or record.get('pending') or record.get('unknown'):
        raise ValueError('A closed linear discovery record is required')
    if any(not np.isfinite(v) or v < 0 for v in rates.values()):
        raise ValueError('Rates must be finite and nonnegative')
    names = list(record['basis'])
    if 'equations_by_channel' in record or 'equations_by_reaction' in record:
        exact = matrix(record)
        required = {str(s) for s in exact.free_symbols}
        allowed = required | {r['rate'] for r in record.get('rates', [])}
        if not required <= rates.keys() or not rates.keys() <= allowed:
            raise ValueError(f'Supply matrix parameters {sorted(required)}')
        values = np.array(exact.subs({sp.Symbol(k): v for k, v in rates.items()}), dtype=float)
    else:
        steps = record['steps']
        channels = {row['channel'] for row in steps}
        if set(rates) != channels:
            raise ValueError(f'Supply channel rates {sorted(channels)}')
        index = {n: i for i, n in enumerate(names)}
        values = np.zeros((len(names), len(names)))
        seen = set()
        for row in steps:
            key = row['source'], row['channel']
            if key in seen:
                raise ValueError('Duplicate source/channel expansion')
            seen.add(key)
            for target, coefficient in row['coordinates'].items():
                values[index[row['source']], index[target]] += rates[row['channel']] * float(Fraction(coefficient))
        if seen != {(name, channel) for name in names for channel in channels}:
            raise ValueError('Incomplete coefficient rows')
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite coefficients')
    return names, values


def trajectory(record, initial, rates, *, end=1., points=201, times=None):
    """Numerical solution of the generated exact expectation system."""
    names, coefficients = linear_matrix(record, rates)
    times = np.asarray(np.linspace(0, end, points) if times is None else times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)):
        raise ValueError('Use increasing finite times starting at zero')
    if set(initial) != set(names):
        raise ValueError('Supply one initial value per discovered coordinate')
    values = np.array([initial[n] for n in names], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Initial values must be finite')
    solution = solve_ivp(lambda t, y: coefficients @ y, (0, times[-1]), values,
                         t_eval=times, method='DOP853', rtol=1e-10, atol=1e-12)
    if not solution.success:
        raise RuntimeError(solution.message)
    return {'time': times, 'values': solution.y.T, 'names': names}
