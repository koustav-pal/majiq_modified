"""
PsiControlsSummary.py

Summary of independent experiments/prefixes quantifications from PsiCoverage
for use as controls

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, Final, List, Sequence, Tuple, Union, cast

import numpy as np
import xarray as xr
from dask.delayed import Delayed
from moccasin.moccasin import persist_with_progress

import rna_majiq.constants as constants

from ..experiments import min_experiments
from .MixinHasEvents import MixinHasEvents
from .PsiCoverage import PsiCoverage
from .PsiGroup import PsiGroup


def _q_from_alpha(alpha: xr.DataArray) -> xr.DataArray:
    """Two-sided quantiles from outer product of alpha, is_lb

    Two-sided quantiles from outer product of alpha, is_lb.
    When is_lb == True, q = 0.5 * alpha. Otherwise, q = 1 - 0.5 * alpha.
    """
    is_lb = xr.DataArray(_ := [False, True], [("is_lb", _)], name="is_lb")
    q_lb = 0.5 * alpha  # lower bound quantiles
    q_ub = 1 - q_lb  # upper bound quantiles
    return q_lb.where(is_lb, q_ub).rename("q")


def _psirange_from_psiquantiles(quantiles: xr.DataArray) -> xr.DataArray:
    """Range between lower/upper quantile (quantiles[..., is_lb])"""
    return (quantiles.sel(is_lb=False) - quantiles.sel(is_lb=True)).rename("psi_range")


class PsiControlsSummary(MixinHasEvents):
    """Summary of PSI posterior means over large group of controls

    Summary of PSI posterior means over large group of controls noting the
    number of experiments which passed and quantiles over the experiments which
    passed

    Parameters
    ----------
    df: xr.Dataset
        With coordinates/variables in `cls.EXPECTED_VARIABLES` and attribute
        `prefixes: List[str]`
    events: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events
    hold_temporary: Sequence[TemporaryDirectory]
        Hold onto temporary directories, preventing their deletion (used by
        :meth:`PsiControlsSummary.from_psi` for temporary store of
        rechunked parameters)

    See Also
    --------
    PsiControlsSummary.from_psi
    PsiControlsSummary.from_zarr
    """

    EXPECTED_VARIABLES: Dict[str, Tuple[str, ...]] = {
        "controls_alpha": ("controls_alpha",),
        "is_lb": ("is_lb",),
        "num_passed": ("ec_idx",),
        "psi_median": ("ec_idx",),
        "psi_quantile": ("ec_idx", "controls_alpha", "is_lb"),
        "total_median": ("ec_idx",),
        "total_quantile": ("ec_idx", "controls_alpha", "is_lb"),
    }

    def __init__(
        self,
        df: xr.Dataset,
        events: xr.Dataset,
        hold_temporary: Sequence[TemporaryDirectory] = tuple(),
    ):
        # save events
        MixinHasEvents.__init__(self, df, events)
        # save df after checking validity
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        if "prefixes" not in df.attrs or not isinstance(df.attrs["prefixes"], list):
            raise ValueError("df.attrs missing required list attribute 'prefixes'")
        # add extra coordinate describing quantiles if not present
        if "q" not in df.variables:
            df = df.assign_coords(controls_q=_q_from_alpha(df["controls_alpha"]))
        self.df: Final[xr.Dataset] = df
        self.hold_temporary: Final = hold_temporary
        return

    @property
    def controls_alpha(self) -> xr.DataArray:
        return self.df["controls_alpha"]

    @property
    def alpha(self) -> xr.DataArray:
        """alias for controls_alpha"""
        return self.controls_alpha

    @property
    def is_lb(self) -> xr.DataArray:
        return self.df["is_lb"]

    @property
    def controls_q(self) -> xr.DataArray:
        return self.df["controls_q"]

    @property
    def q(self) -> xr.DataArray:
        """alias for controls_q"""
        return self.controls_q

    @property
    def num_passed(self) -> xr.DataArray:
        return self.df["num_passed"]

    @property
    def prefixes(self) -> List[str]:
        return self.df.attrs["prefixes"]

    @property
    def num_prefixes(self) -> int:
        return len(self.prefixes)

    def passed_min_experiments(
        self,
        min_experiments_f: Union[
            float, Sequence[float]
        ] = constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
    ) -> xr.DataArray:
        """Get boolean mask of events that pass enough experiments"""
        if isinstance(min_experiments_f, float):
            min_experiments_f = [min_experiments_f]
        self_min_experiments = min_experiments(
            xr.DataArray(min_experiments_f, [("min_experiments_f", min_experiments_f)]),
            self.num_prefixes,
        )
        return self.num_passed >= self_min_experiments

    @property
    def psi_median(self) -> xr.DataArray:
        return self.df["psi_median"]

    @property
    def psi_quantile(self) -> xr.DataArray:
        return self.df["psi_quantile"]

    @property
    def total_median(self) -> xr.DataArray:
        return self.df["total_median"]

    @property
    def total_quantile(self) -> xr.DataArray:
        return self.df["total_quantile"]

    @cached_property
    def psi_range(self) -> xr.DataArray:
        """For each controls_alpha, range between lower/upper quantiles (scale)"""
        return _psirange_from_psiquantiles(self.psi_quantile)

    @classmethod
    def from_psi(
        cls,
        psi: Union[PsiGroup, PsiCoverage],
        alpha: Union[float, Sequence[float]] = constants.DEFAULT_OUTLIERS_ALPHA,
    ) -> "PsiControlsSummary":
        if isinstance(alpha, float):
            alpha = [alpha]
        else:
            # make sure alpha is sorted and has unique elements
            alpha = sorted(set(alpha))
        # get everything for df but the quantiles
        df = xr.Dataset(
            dict(),
            dict(
                num_passed=psi.num_passed,
                controls_alpha=("controls_alpha", alpha),
            ),
            dict(prefixes=psi.prefixes),
        ).assign_coords(
            controls_q=lambda df: _q_from_alpha(df["controls_alpha"]),
        )
        stacked_q = df["controls_q"].stack(_idx=["controls_alpha", "is_lb"])
        # add median to quantiles
        compute_q = xr.concat(
            [
                xr.DataArray(np.array([0.5], dtype=stacked_q.dtype), dims=["_idx"]),
                stacked_q.drop_vars(["_idx", "is_lb", "controls_alpha", "controls_q"]),
            ],
            dim="_idx",
        )
        computed_quantiles = psi.raw_psi_mean_population_quantile(
            compute_q.values.tolist(), "controls_q"
        ).swap_dims(controls_q="_idx")
        computed_totals = psi.raw_total_population_quantile(
            compute_q.values.tolist(), "controls_q"
        ).swap_dims(controls_q="_idx")
        psi_median = computed_quantiles.isel(_idx=0, drop=True)
        psi_quantile = (
            computed_quantiles.isel(_idx=slice(1, None))
            .assign_coords(_idx=stacked_q["_idx"])
            .unstack("_idx")
        )
        total_median = computed_totals.isel(_idx=0, drop=True)
        total_quantile = (
            computed_totals.isel(_idx=slice(1, None))
            .assign_coords(_idx=stacked_q["_idx"])
            .unstack("_idx")
        )
        df = df.assign(
            psi_median=psi_median,
            psi_quantile=psi_quantile,
            total_median=total_median,
            total_quantile=total_quantile,
        )
        hold_temporary = psi.hold_temporary if isinstance(psi, PsiGroup) else tuple()
        return PsiControlsSummary(
            df.chunk({"ec_idx": -1}),
            psi.events_df,
            hold_temporary=hold_temporary,
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "PsiControlsSummary":
        return PsiControlsSummary(
            xr.open_zarr(path, group=constants.NC_PSICONTROLS),
            xr.open_zarr(path, group=constants.NC_EVENTS),
        )

    def to_zarr(
        self,
        path: Union[str, Path],
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save PSI coverage dataset as zarr

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        consolidated: bool
            When saving the file make sure that it is consolidated. In general,
            if you are appending a bunch of files together, it can make sense
            to set consolidated=False, and consolidate on the last write (only
            consolidate once). But, don't forget to consolidate at the end.
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        # save df
        save_df_future = cast(
            Delayed,
            self.df.chunk(self.df.sizes).to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICONTROLS,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = persist_with_progress(save_df_future)
        save_df_future.compute()
        # save events
        self.events_to_zarr(path, mode="a", consolidated=consolidated)
        return

    def __repr__(self) -> str:
        MAX_PREFIXES_END = 1  # how many prefixes on either end to display
        if self.num_prefixes > 2 * MAX_PREFIXES_END:
            print_prefixes_list = [
                *self.prefixes[:MAX_PREFIXES_END],
                *([] if self.num_prefixes <= 2 * MAX_PREFIXES_END else ["..."]),
                *self.prefixes[-MAX_PREFIXES_END:],
            ]
        else:
            print_prefixes_list = self.prefixes
        print_prefixes = ", ".join(print_prefixes_list)
        return (
            f"PsiControlsSummary[{self.num_connections}]"
            f" for {self.num_prefixes} experiments [{print_prefixes}]"
        )
