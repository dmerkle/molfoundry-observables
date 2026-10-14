"""The public benchmark accepts zero warm-ups and checks every timed result."""
import pytest

from scripts.comparison.ab_scaling import molfoundry_benchmark as benchmark


def test_zero_warmups_counts_the_same_compositions():
    cold = benchmark.benchmark('alt', 2, warmups=0, repetitions=2)
    warm = benchmark.benchmark('alt', 2, warmups=1, repetitions=1)
    assert cold['raw_composites'] == warm['raw_composites'] == 21
    assert cold['per_order'] == warm['per_order']
    assert len(cold['times_seconds']) == 2


def test_zero_warmups_still_rejects_changed_results(monkeypatch):
    results = iter([{'raw_composites': 21}, {'raw_composites': 22}])
    monkeypatch.setattr(benchmark, 'enumerate_raw', lambda rows: next(results))
    with pytest.raises(RuntimeError, match='non-deterministic'):
        benchmark.benchmark('alt', 2, warmups=0, repetitions=2)
