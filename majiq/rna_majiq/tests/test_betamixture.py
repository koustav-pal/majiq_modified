"""
test_betamixture.py

Test rna_majiq beta_mixture implementations

Author: Joseph K Aicher
"""

import rna_majiq.beta_mixture as bm
import numpy as np
import pytest


@pytest.mark.parametrize("q", [0.01, 0.1, 0.5, 0.9, 0.95])
@pytest.mark.parametrize("psi", [0.1, 0.24, 0.26, 0.9])
@pytest.mark.parametrize("total", [1, 100, 10000])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_betamixture_quantile(q, psi, total, dtype):
    # XXX: we currently have a bug which causes overflow warning when
    #  1 < a < b/4 (or 1 < b < a/4)
    a = psi * total
    b = total - a
    # this should evaluate without any warnings
    bm.quantile(
        np.array(q, dtype=dtype), np.array([a], dtype=dtype), np.array([b], dtype=dtype)
    )
    return
