"""
FactorsModel.py

Wrapper of moccasin.ModelUnknownConfounders

Author: Joseph K Aicher
"""

import logging
from pathlib import Path
from typing import Final, Optional, Union

import moccasin.moccasin as mc
import numpy as np
import xarray as xr
from moccasin.moccasin import persist_with_progress

import rna_majiq.constants as constants

from ..core.PsiCoverage import PsiCoverage
from ..logger import get_logger
from .Factors import Factors


class FactorsModel(object):
    """Wrapper of `moccasin.ModelUnknownConfounders` for :class:`PsiCoverage`"""

    def __init__(self, model: mc.ModelUnknownConfounders):
        self._model: Final[mc.ModelUnknownConfounders] = model
        return

    def predict_unknown(
        self,
        psicov: PsiCoverage,
        factors: Factors,
        num_unknown_factors: Optional[int] = None,
        show_progress: bool = False,
    ) -> Factors:
        """Create :class:`Factors` for unknown factors using :class:`PsiCoverage`"""
        try:
            factors = factors[psicov.prefixes]
        except KeyError as err:
            missing = set(psicov.prefixes) - set(factors.prefixes)
            raise ValueError(f"Missing input factors for prefixes: {missing}") from err
        factors_arr = self.model.predict(
            psicov.raw_psi,
            psicov.event_passed,
            factors.factors,
            num_factors=num_unknown_factors,
        )
        if show_progress:
            factors_arr = persist_with_progress(factors_arr)
        factors_arr.load()
        return Factors(
            factors_arr.rename("factors")
            .to_dataset()
            .rename_dims(new_factor="factor")
            .rename_vars(new_factor="factor")
        )

    def predict(
        self,
        psicov: PsiCoverage,
        factors: Factors,
        num_unknown_factors: Optional[int] = None,
        show_progress: bool = False,
    ) -> Factors:
        """Create updated :class:`Factors` using :class:`PsiCoverage`"""
        return factors.combine(
            self.predict_unknown(
                psicov,
                factors,
                num_unknown_factors=num_unknown_factors,
                show_progress=show_progress,
            )
        )

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Save :class:`FactorsModel` to specified path"""
        return self.model.to_zarr(path)

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "FactorsModel":
        return FactorsModel(mc.ModelUnknownConfounders.from_zarr(path))

    @classmethod
    def train(
        cls,
        psicov: PsiCoverage,
        factors: Factors,
        default_num_factors: int = constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
        max_events: int = constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        log: Optional[logging.Logger] = None,
        show_progress: bool = False,
    ) -> "FactorsModel":
        if log is None:
            log = get_logger()
        try:
            factors = factors[psicov.prefixes]
        except KeyError as err:
            missing = set(psicov.prefixes) - set(factors.prefixes)
            raise ValueError(f"Missing input factors for prefixes: {missing}") from err
        model = mc.ModelUnknownConfounders.train(
            psicov.raw_psi,
            psicov.event_passed,
            psicov.lsv_offsets.values.view(np.int64),
            factors.factors,
            default_num_factors=default_num_factors,
            max_events=max_events,
            dim_prefix="prefix",
            dim_factor="factor",
            dim_ecidx="ec_idx",
            log=log,
            show_progress=show_progress,
        )
        return FactorsModel(model)

    @property
    def model(self) -> mc.ModelUnknownConfounders:
        """Underlying Moccasin model of unknown confounders"""
        return self._model

    @property
    def explained_variance(self) -> xr.DataArray:
        """proportion of training data residuals variance explained by each new factor"""
        return self.model.explained_variance

    def explained_variance_report(self) -> str:
        """String explanation of variance explained by model and per factor"""
        explained_variance = self.explained_variance.values
        return (
            f"Learned {self.max_num_factors} factors explaining"
            f" {explained_variance.sum():.3%} of residual variance"
            " (after accounting for known factors):"
            + "".join(
                [
                    f"\nfactor {i}\t{x:.3%} (cumulative: {cum_x:.3%})"
                    for i, (x, cum_x) in enumerate(
                        zip(explained_variance, explained_variance.cumsum()), 1
                    )
                ]
            )
        )

    @property
    def max_num_factors(self) -> int:
        """Maximum number of new factors to infer limited by stored components"""
        return self.model.max_num_factors
