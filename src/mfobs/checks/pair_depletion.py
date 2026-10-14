"""Independent, exact generator audit for a balanced pair-consuming channel.

This is not a closure reducer. The polynomial certificate and finite batch
matrices are verification outputs, never input to molecular composition.
"""
from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix

N = sp.Symbol('n', integer=True, nonnegative=True)
K = sp.Symbol('k', positive=True)


def generator(polynomial):
    """n -> n-2 at total propensity k*binomial(n, 2)."""
    f = sp.sympify(polynomial)
    return sp.expand(K*N*(N-1)/2 * (f.subs(N, N-2) - f))


def polynomial_certificate(order=8):
    """Check finitely many instances of the separately documented induction.

    Finite checks alone are NOT the all-order proof. The induction is that
    Q(c*n**d + lower) = -k*d*c*n**(d+1) + lower, for d >= 1.
    """
    if type(order) is not int or order < 1:
        raise ValueError('order must be a positive integer')
    f, rows = N, []
    for r in range(order+1):
        poly = sp.Poly(f, N)
        leading = sp.expand((-K)**r * sp.factorial(r))
        assert poly.degree() == r+1
        assert sp.expand(poly.LC()-leading) == 0
        rows.append({'iteration': r, 'polynomial': str(f),
                     'degree': poly.degree(), 'leading_coefficient': str(leading)})
        f = generator(f)
    return rows


def finite_batch(n0):
    """Exact rational Krylov rank for one finite batch, using unit rate.

    States contain only unreacted ethylene and inert cyclobutane. The rank
    is for the observable n, not the rank of the entire generator matrix.
    A separately supplied constant coordinate is not counted.
    """
    if type(n0) is not int or n0 < 0:
        raise ValueError('n0 must be a nonnegative integer')
    states = list(range(n0, -1, -2))
    q = sp.zeros(len(states))
    for i, n in enumerate(states[:-1]):
        rate = sp.Rational(n*(n-1), 2)
        q[i, i], q[i, i+1] = -rate, rate
    v, columns = sp.Matrix(states), []
    for _ in states:
        columns.append(v)
        v = q*v
    # Use exact domain arithmetic: generic expression elimination grows huge
    # intermediate integers for the larger batch audits.
    rank = DomainMatrix.from_Matrix(sp.Matrix.hstack(*columns)).rank()
    return {'initial_ethylene': n0, 'reachable_population_states': len(states),
            'states': states, 'seed_krylov_dimension': rank,
            'rate': 'k=1 (positive scalar does not change rank)'}
