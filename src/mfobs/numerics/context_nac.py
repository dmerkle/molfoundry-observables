"""Finite forward equations of the generated chlorination event model.

The reachable coordinate states are explored only after symbolic synthesis.
No ring masks or graph transition table are used to build this generator.
"""
from collections import deque
import numpy as np
import sympy as sp
from scipy.linalg import expm

from mfobs.numerics.events import polynomial


def initial_counts(result, mask=0):
    from mfobs.core.pattern_support import Pattern
    from mfobs.core.engine import mf
    from mfobs.models.chlorination.inputs import query_patterns
    patterns = query_patterns()
    for name, record in result['synthesis']['positive_patterns'].items():
        if name not in patterns:
            patterns[name] = Pattern(name, mf.Rule.fromGMLString(record['gml'], add=False), 'positive')
    return _coordinate_values(mask, patterns)


def forward_system(result, initial, maximum_states=100):
    if result['status'] != 'PASS' or result['pjc']['status'] != 'CLOSED_FINITE_EVENT_MODEL':
        raise ValueError('A certified event model is required')
    names = tuple(result['pjc']['coordinates'])
    if set(initial) != set(names) or any(type(v) is not int or v < 0 for v in initial.values()):
        raise ValueError('Supply an integer initial count for every coordinate')
    symbols = {n: sp.Symbol(n) for n in names}
    activities = [polynomial(e['activity'], symbols) for e in result['events']]
    deltas = [tuple(e['delta'].get(n, 0) for n in names) for e in result['events']]
    start = tuple(initial[n] for n in names)
    states, queue, transitions = {start: 0}, deque([start]), []
    while queue:
        state = queue.popleft()
        for activity, delta in zip(activities, deltas):
            rate = activity.subs(dict(zip(symbols.values(), state)))
            if not rate.is_Rational or rate < 0:
                raise ValueError(f'Invalid activity {rate} at {state}')
            if not rate:
                continue
            target = tuple(x + d for x, d in zip(state, delta))
            if any(x < 0 for x in target):
                raise ValueError('An enabled jump produces negative counts')
            if target not in states:
                if len(states) >= maximum_states:
                    raise RuntimeError('Reachability budget exhausted; no truncated solution returned')
                states[target] = len(states)
                queue.append(target)
            transitions.append((states[state], states[target], float(rate)))
    generator = np.zeros((len(states), len(states)))
    for source, target, rate in transitions:
        generator[source, target] += rate
        generator[source, source] -= rate
    return names, tuple(states), generator


def trajectory(result, initial, *, times, rate=1.):
    times = np.asarray(times, float)
    if (not np.isfinite(rate) or rate < 0 or times.ndim != 1 or not len(times)
            or times[0] != 0 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)):
        raise ValueError('Use a finite nonnegative rate and increasing times starting at zero')
    names, states, generator = forward_system(result, initial)
    probabilities = np.asarray([expm(rate*t*generator)[0] for t in times])
    if np.min(probabilities) < -1e-10 or not np.allclose(probabilities.sum(1), 1., atol=1e-10):
        raise ArithmeticError('Forward equation probability check failed')
    z = np.asarray([s[names.index('Z')] for s in states], float)
    mean = probabilities @ z
    variance = probabilities @ (z*z) - mean*mean
    if np.min(variance) < -1e-10:
        raise ArithmeticError('Negative variance')
    incomplete = z < result['synthesis']['state_language']['chlorine_total']
    blocked = (np.diag(generator) == 0) & incomplete
    return {'time': times, 'Z': mean, 'variance': np.maximum(variance, 0.),
            'blocked': probabilities @ blocked, 'states': states,
            'generator': generator, 'probabilities': probabilities}




def _coordinate_values(mask, patterns, *, n=0, methane=0, h2=0, hcl=0):
    """Count the generated positive observable embeddings with MolFoundry."""
    from mfobs.core.molfoundry_execution import count_matches
    from mfobs.models.chlorination.inputs import host as molecular_host
    host = molecular_host(n=n, mask=mask, methane=methane, hcl=hcl, h2=h2)
    return {name: count_matches(pattern.rule, host)
            for name, pattern in patterns.items()}
