import pytest
import sympy as sp

from mfobs.checks.pair_depletion import N, K, generator, polynomial_certificate, finite_batch


def test_first_two_generators_without_factorizing_expectations():
    assert generator(N) == -K*N**2 + K*N
    assert sp.expand(generator(N**2) - (-2*K*N**3+4*K*N**2-2*K*N)) == 0


def test_degree_growth_and_exact_leading_coefficients():
    rows = polynomial_certificate(8)
    assert [r['degree'] for r in rows] == list(range(1, 10))


@pytest.mark.parametrize('n0', [0, 1, 2, 3, 4, 5, 8, 9, 16, 17, 32, 33])
def test_fixed_batch_rank_grows_with_inventory(n0):
    result = finite_batch(n0)
    expected = (n0+1)//2
    assert result['seed_krylov_dimension'] == expected
    assert result['reachable_population_states'] == n0//2+1


@pytest.mark.parametrize('n0', [-1, True, 1.5])
def test_invalid_batch_size(n0):
    with pytest.raises(ValueError):
        finite_batch(n0)
