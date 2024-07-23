"""
test_stats.py

Test rna_majiq statistics implementations

Author: Joseph K Aicher
"""

import itertools

import rna_majiq.stats as stats
import numpy as np
import pytest
from scipy.special import comb


@pytest.fixture
def random_state() -> np.random.RandomState:
    SEED: int = 11402022
    return np.random.RandomState(SEED)


@pytest.mark.parametrize("test", ["mannwhitneyu", "tnom", "infoscore"])
@pytest.mark.parametrize("n", [2, 3, 4, 10, 50, 1000])
def test_stats_nonparametric_identical(
    random_state: np.random.RandomState, test: str, n: int
):
    assert n > 1, "invalid test"
    x = np.ones(n)
    sortx = np.arange(n)
    # do we enumerate all possible labels or just sample
    MAX_REPEATS = 1024
    if (1 << n) <= MAX_REPEATS:
        labels = np.array([list(y) for y in itertools.product([True, False], repeat=n)])
    else:
        labels = random_state.choice([True, False], size=(MAX_REPEATS, n))
    # ignore labels with all True or all False
    labels = labels[labels.any(axis=-1) & (~labels.any(axis=-1))]
    # compute pvalues
    pvalues = getattr(stats, test)(x, labels, sortx=sortx)
    # these should all be 1
    np.testing.assert_equal(pvalues, 1)
    return


@pytest.mark.parametrize("test", ["mannwhitneyu", "tnom", "infoscore"])
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 10, 50, 1000])
def test_stats_nonparametric_separated(test: str, n: int):
    """permutation tests when groups are linearly separable closed form"""
    if test == "mannwhitneyu" and n > 64:
        pytest.skip(reason="normal approximation with tie-correction != permutation")
        return
    if test == "infoscore" and n > 100:
        pytest.skip(reason="infoscore with large n becomes too slow")
        return
    # all possible n1/n2 (including n1 = 0 or n2 = 0)
    labels = np.triu(np.ones((1 + n, n), dtype=bool))
    x = np.where(labels, 1.0, 0.0)
    sortx = np.arange(n)
    # compute pvalues
    pvalues = getattr(stats, test)(x, labels, sortx=sortx)
    # expected value?
    n1 = labels.sum(axis=-1)
    n2 = (~labels).sum(axis=-1)
    expected = np.where((n1 > 0) & (n2 > 0), 2 / comb(n, n1), np.nan)
    np.testing.assert_almost_equal(pvalues, expected)
    return


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 10, 50, 1000])
@pytest.mark.parametrize("identical", [False, True])
def test_stats_ttest_zero_variance(n: int, identical: bool):
    # all possible n1/n2 (including n1 = 0 or n2 = 0)
    labels = np.triu(np.ones((1 + n, n), dtype=bool))
    x = np.where(identical | labels, 1.0, 0.0)
    # compute pvalues
    pvalues = stats.ttest(x, labels)
    # expected value?
    n1 = labels.sum(axis=-1)
    n2 = (~labels).sum(axis=-1)
    expected = np.where((n1 > 1) & (n2 > 1), 1 if identical else 0, np.nan)
    np.testing.assert_equal(pvalues, expected)
    return


@pytest.mark.parametrize("test", ["ttest", "mannwhitneyu", "tnom", "infoscore"])
@pytest.mark.parametrize("n", [1, 2, 3, 4, 10, 50, 1000])
def test_stats_single_label(random_state: np.random.RandomState, test: str, n: int):
    REPEATS = 1024
    x = random_state.normal(size=(REPEATS, 1, n))
    labels = np.array([[True] * n, [False] * n])
    pvalues = getattr(stats, test)(x, labels)
    # these should all be nan
    np.testing.assert_equal(pvalues, np.nan)
    return
