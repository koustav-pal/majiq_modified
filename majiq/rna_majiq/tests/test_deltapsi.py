"""
test_deltapsi.py

Test rna_majiq deltapsi implementations

Author: Joseph K Aicher
"""

import rna_majiq as nm
import numpy as np
import pytest


def test_deltapsi_prior_update_zero():
    """tests that DPsiPrior can do update with zero events"""
    psi1 = nm.PsiCoverage.mock_with_psi_and_total([], [], "grp1")
    psi2 = nm.PsiCoverage.mock_with_psi_and_total([], [], "grp2")
    init_prior = nm.DPsiPrior()
    update_prior = init_prior.empirical_update(psi1, psi2)
    assert init_prior == update_prior
    return


@pytest.mark.parametrize("n_update_a", [0, 1, 3])
@pytest.mark.parametrize("legacy", [False, True])
def test_deltapsi_prior_update_thresholds(n_update_a, legacy):
    """tests that DPsiPrior does/does not update given thresholds"""
    minreads = 30.0
    nlsvs = 100
    coverage = (1 + minreads) * np.power(10, np.linspace(0, 6, nlsvs))
    psi1 = np.random.beta(0.5, 0.5, size=nlsvs).astype(np.float32)
    psi1 = np.where(psi1 < 0.5, 1 - psi1, psi1)
    total1 = coverage / psi1
    # NOTE: completely independent values of psi2 vs psi1
    # this is not testing what the resulting priors look like
    psi2 = np.random.beta(0.5, 0.5, size=nlsvs).astype(np.float32)
    psi2 = np.where(psi2 < 0.5, 1 - psi2, psi2)
    total2 = coverage / psi2
    init_prior = nm.DPsiPrior()
    psicov1 = nm.PsiCoverage.mock_with_psi_and_total(psi1, total1, "grp1")
    psicov2 = nm.PsiCoverage.mock_with_psi_and_total(psi2, total2, "grp2")
    assert init_prior != init_prior.empirical_update(
        psicov1,
        psicov2,
        minreads=minreads,
        min_lsvs=nlsvs,
        n_update_a=n_update_a,
        legacy=legacy,
    )
    assert init_prior == init_prior.empirical_update(
        psicov1,
        psicov2,
        minreads=2 + minreads,
        min_lsvs=nlsvs,
        n_update_a=n_update_a,
        legacy=legacy,
    )
    assert init_prior == init_prior.empirical_update(
        psicov1,
        psicov2,
        minreads=minreads,
        min_lsvs=nlsvs + 1,
        n_update_a=n_update_a,
        legacy=legacy,
    )
    return


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("eps_max", [0.0, 1e-5, 1e-2, 1e-1])
def test_deltapsi_prior_update_nonsingular(legacy, eps_max):
    """When values of psi are identical, can get singularity in prior update"""
    minreads = 30.0
    nlsvs = 100
    coverage = minreads * np.power(10, np.linspace(0, 6, nlsvs))
    psi1 = np.random.beta(0.5, 0.5, size=nlsvs).astype(np.float32)
    psi1 = np.where(psi1 < 0.5, 1 - psi1, psi1)
    total1 = coverage / psi1
    psi2 = psi1 + np.random.uniform(0, np.minimum(1 - psi1, eps_max))
    init_prior = nm.DPsiPrior()
    psicov1 = nm.PsiCoverage.mock_with_psi_and_total(psi1, total1, "grp1")
    psicov2 = nm.PsiCoverage.mock_with_psi_and_total(psi2, total1, "grp2")
    updated_prior = init_prior.empirical_update(
        psicov1, psicov2, minreads=minreads, min_lsvs=nlsvs, legacy=legacy
    )
    assert (
        np.isfinite(updated_prior.pmix).all() and np.isfinite(updated_prior.a).all()
    ), f"{updated_prior} values must all be finite"
    np.testing.assert_almost_equal(
        updated_prior.pmix.sum(),
        1.0,
        err_msg=f"{updated_prior} must have pmix summing to 1",
    )
    return
