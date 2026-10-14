"""Analytical probability laws and explicit failure cases for the FSP consumer."""
import numpy as np
import pytest
import sympy as sp

from mfobs.numerics.finite_state import solve_moments


def immigration_consumption(*, influx=2., initial_count=3, tolerance=1e-6,
                            duration=3., death=.7, adaptive=None):
    # A different model from the word example: each token is consumed once,
    # increasing a passive counter. Its exact law is Binomial + Poisson.
    times = np.linspace(0., duration, 31)
    x = sp.Symbol('available')
    solution = solve_moments(
        names=('count', 'available'), deltas=((0, 1), (1, -1)),
        propensities=(sp.Float(influx), death*x),
        initial={'count': 2, 'available': initial_count}, query='count',
        boolean_names=(), times=times, truncation_tolerance=tolerance, adaptive=adaptive,
    )
    success = 1-np.exp(-death*times)
    poisson_mean = influx*(times-success/death)
    mean = 2+initial_count*success+poisson_mean
    variance = initial_count*success*(1-success)+poisson_mean
    return solution, mean, variance


def test_finite_consumption_matches_binomial_law_without_truncation():
    solution, mean, variance = immigration_consumption(influx=0.)
    np.testing.assert_allclose(solution.mean, mean, atol=2e-12)
    np.testing.assert_allclose(solution.variance, variance, atol=2e-12)
    np.testing.assert_allclose(solution.surviving_probability, 1., atol=2e-12)
    assert solution.diagnostics['exit_probability_bound'] == 0
    assert solution.diagnostics['states'] == 4


@pytest.mark.parametrize('tolerance', [1e-6, .25])
def test_infinite_model_moment_bounds_cover_known_law(tolerance):
    solution, mean, variance = immigration_consumption(tolerance=tolerance)
    info = solution.diagnostics
    assert np.max(np.abs(solution.mean-mean)) <= info['mean_error_bound'] + 1e-11
    assert np.max(np.abs(solution.variance-variance)) <= info['variance_error_bound'] + 1e-10
    assert info['variance_error_bound'] <= tolerance
    assert info['observed_mass_defect'] <= info['exit_probability_bound'] + 1e-12
    # The envelope is discovered from the event jumps, without pool/query names.
    assert info['query_envelope']['weights'] == {'count': 1, 'available': 1}


def test_pure_poisson_count_has_no_base_state_and_handles_self_transitions():
    times = np.linspace(0, 2, 11)
    solution = solve_moments(names=('N',), deltas=((1,),), propensities=(sp.Integer(3),),
                             initial={'N': 4}, query='N', boolean_names=(), times=times)
    np.testing.assert_allclose(solution.mean, 4+3*times, atol=1e-12)
    np.testing.assert_allclose(solution.variance, 3*times, atol=1e-12)
    assert solution.diagnostics['states'] == 1


def test_zero_rates_leave_query_deterministic():
    s, mean, variance = immigration_consumption(influx=0., initial_count=0)
    np.testing.assert_allclose(s.mean, mean, atol=1e-12)
    np.testing.assert_allclose(s.variance, variance, atol=1e-12)


def test_probability_is_killed_at_boundary_not_reflected_or_renormalised():
    from scipy.stats import poisson
    solution = solve_moments(
        names=('q', 'x'), deltas=((1, 1),), propensities=(sp.Integer(2),),
        initial={'q': 0, 'x': 0}, query='q', boolean_names=(),
        times=np.linspace(0, 3, 31), truncation_tolerance=500.,
    )
    cutoff = solution.diagnostics['cutoff']['x']
    counts = np.arange(cutoff+1)
    probabilities = poisson.pmf(counts, 6.)
    np.testing.assert_allclose(solution.surviving_probability[-1], probabilities.sum(), atol=1e-12)
    np.testing.assert_allclose(solution.mean[-1], counts @ probabilities, atol=1e-12)
    np.testing.assert_allclose(solution.second_moment[-1], counts**2 @ probabilities, atol=1e-12)
    np.testing.assert_allclose(solution.diagnostics['mean_error_bound'], 6-counts @ probabilities, atol=1e-11)
    np.testing.assert_allclose(solution.diagnostics['second_moment_error_bound'], 42-counts**2 @ probabilities, atol=1e-10)
    assert solution.surviving_probability[-1] < .9
    assert np.all(np.diff(solution.surviving_probability) <= 1e-12)


def test_adaptive_projection_grows_and_bounds_both_moments_at_long_time():
    solution, mean, variance = immigration_consumption(duration=100., death=.01, adaptive=True)
    info = solution.diagnostics
    assert len(info['projection_attempts']) > 1
    assert info['variance_error_bound'] <= 1e-6
    assert np.max(np.abs(solution.mean-mean)) <= info['mean_error_bound'] + 1e-9
    assert np.max(np.abs(solution.variance-variance)) <= info['variance_error_bound'] + 1e-8


def test_query_dependent_rates_and_unbounded_branching_are_rejected():
    q, x = sp.symbols('q x')
    kwargs = dict(names=('q', 'x'), deltas=((1, -1), (0, 1)), initial={'q': 0, 'x': 1},
                  query='q', boolean_names=(), times=[0, 1])
    with pytest.raises(ValueError, match='absent from all propensities'):
        solve_moments(**kwargs, propensities=(q*x, sp.Integer(1)))
    with pytest.raises(ValueError, match='state-dependent event increases'):
        solve_moments(**kwargs, propensities=(x, x))
    with pytest.raises(ValueError, match='No query-moment envelope'):
        solve_moments(**{**kwargs, 'deltas': ((1, 0), (0, 1))}, propensities=(x, sp.Integer(1)))


def test_budget_and_invalid_domains_fail_instead_of_returning_a_partial_curve():
    kwargs = dict(names=('q', 'x'), deltas=((1, -1),), propensities=(sp.Symbol('x'),),
                  initial={'q': 0, 'x': 2}, query='q', boolean_names=(), times=[0, 1])
    with pytest.raises(RuntimeError, match='state budget'):
        solve_moments(**kwargs, max_states=1)
    with pytest.raises(ValueError, match='declared coordinate domain'):
        solve_moments(**{**kwargs, 'propensities': (sp.Integer(1),)})
    with pytest.raises(ValueError, match='uniform time grid'):
        solve_moments(**{**kwargs, 'times': [0, .2, 1]})
    with pytest.raises(ValueError, match='integer jump'):
        solve_moments(**{**kwargs, 'deltas': ((.5, -1),)})
