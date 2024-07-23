"""
CoverageModel.py

Wrapper of moccasin functions to model and predict :class:`PsiCoverage`

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Dict, Final, MutableMapping, Optional, Tuple, Union, cast

import dask.array as da
import moccasin.moccasin as mc
import numpy as np
import xarray as xr
from dask.delayed import Delayed
from moccasin.moccasin import persist_with_progress

import rna_majiq.constants as constants
from rna_majiq._offsets import clip_and_normalize_strict

from ..core.PsiCoverage import PsiCoverage
from .Factors import Factors


class CoverageModel(object):
    """Wrapper over Moccasin functionality for :class:`PsiCoverage`

    Parameters
    ----------
    df: xr.Dataset
        With variables as in `CoverageModel.EXPECTED_VARIABLES`
    """

    EXPECTED_VARIABLES: Final[Dict[str, Tuple[str, ...]]] = {
        "factor": ("factor",),
        "raw_model": (
            "ec_idx",
            "factor",
        ),
        "bootstrap_model": ("ec_idx", "bootstrap_replicate", "factor"),
    }
    # ZARR_GROUP = "moccasin_coverage_model"
    ZARR_GROUP = None

    def __init__(self, df: xr.Dataset):
        # save df after checking validity
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        # clear any encodings on variables (i.e. original zarr chunks)
        for v in df.variables.values():
            v.encoding.clear()
        # do not load into memory
        self.df: Final[xr.Dataset] = df
        return

    @property
    def raw_model(self) -> xr.DataArray:
        """Model parameters for raw coverage"""
        return self.df["raw_model"]

    @property
    def bootstrap_model(self) -> xr.DataArray:
        """Model parameters for bootstrap coverage"""
        return self.df["bootstrap_model"]

    def to_zarr(
        self, store: Union[MutableMapping, str, Path], show_progress: bool = False
    ) -> None:
        """Save :class:`CoverageModel` as zarr to specified path/store.

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            Path for output file/store in zarr format
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        # chunk over ec_idx
        use_chunks = dict(self.df.sizes)
        use_chunks["ec_idx"] = constants.DEFAULT_COVERAGE_CHUNKS
        save_df = self.df.chunk(use_chunks)
        # save result
        save_df_future = cast(
            Delayed,
            save_df.to_zarr(store, mode="w", group=self.ZARR_GROUP, compute=False),
        )
        if show_progress:
            save_df_future = persist_with_progress(save_df_future)
        save_df_future.compute()
        return

    @classmethod
    def from_zarr(cls, store: Union[MutableMapping, str, Path]) -> "CoverageModel":
        """Load :class:`CoverageModel` from zarr at specified path/store"""
        return CoverageModel(xr.open_zarr(store, group=cls.ZARR_GROUP))

    @staticmethod
    def _fused_raw_and_bootstrap(
        raw_arr: xr.DataArray, bootstrap_arr: xr.DataArray
    ) -> xr.DataArray:
        """Fuse raw and bootstrap arrays together for shared computing.

        Fuse raw and bootstrap arrays together for shared computing. Expects
        that arrays have same dimensions except "bootstrap_replicate" in
        bootstrap_arr. Returns array with one more element in bootstrap
        replicate, where raw array is the first value.
        """
        bootstrap_sizes = dict(bootstrap_arr.sizes)
        try:
            bootstrap_sizes.pop("bootstrap_replicate")
        except KeyError as err:
            raise ValueError(
                "bootstrap_arr needs to have bootstrap_replicate dimension"
            ) from err
        if "bootstrap_replicate" in raw_arr.dims:
            raise ValueError("raw_arr must not have bootstrap_replicate dimension")
        if raw_arr.sizes != bootstrap_sizes:
            raise ValueError(
                "raw_arr and bootstrap_arr do not agree on dimension sizes"
                " excluding bootstrap_replicate:"
                f" ({raw_arr.sizes = }, {bootstrap_arr.sizes = })"
            )
        return xr.concat(
            [raw_arr.expand_dims(bootstrap_replicate=1), bootstrap_arr],
            dim="bootstrap_replicate",
        ).chunk({"bootstrap_replicate": -1})

    @staticmethod
    def _raw_from_fused(fused_arr: xr.DataArray) -> xr.DataArray:
        """Extract first bootstrap from fused array (corresponding to raw)"""
        return fused_arr.isel(bootstrap_replicate=0, drop=True)

    @staticmethod
    def _bootstrap_from_fused(fused_arr: xr.DataArray) -> xr.DataArray:
        """Extract all but first bootstrap from fused array (corresponding to bootstrap)"""
        return fused_arr.isel(bootstrap_replicate=slice(1, None))

    @classmethod
    def train(
        self,
        psicov: PsiCoverage,
        factors: Factors,
        compute_chunksize: Optional[int] = None,
    ) -> "CoverageModel":
        """Train :class:`CoverageModel` on raw and bootstrapped coverage

        Parameters
        ----------
        psicov: PsiCoverage
            Prefixes with coverage to model.
            Generally, should be opened with ec_idx_nchunks=1 to avoid memory
            issues.
        factors: Factors
            Confounding/non-confounding factors for modeling coverage
        compute_chunksize: Optional[int]
            Rechunk connections by this chunksize for the computation.
            Prefixes will all be in a single chunk.
            If not specified, MAJIQ will automatically pick a chunksize.
        """
        # subset factors to match psicov prefixes
        try:
            factors = factors[psicov.prefixes]
        except KeyError as err:
            missing = set(psicov.prefixes) - set(factors.prefixes)
            raise ValueError(f"Missing input factors for prefixes: {missing}") from err
        if psicov.num_events > 0:
            # solve model for raw/bootstrap psi simultaneously
            compute_psi = self._fused_raw_and_bootstrap(
                psicov.raw_psi, psicov.bootstrap_psi
            )
            compute_model = mc.infer_model_params(
                compute_psi,
                psicov.event_passed,
                factors.factors,
                complete=False,
                dim_prefix="prefix",
                dim_factor="factor",
            )
        else:
            compute_model = xr.DataArray(
                da.empty(
                    (0, factors.num_factors, 1 + psicov.num_bootstraps),
                    dtype=psicov.bootstrap_psi.dtype,
                ),
                coords={"factor": factors.factor_names},
                dims=["ec_idx", "factor", "bootstrap_replicate"],
            )
        raw_model = self._raw_from_fused(compute_model)
        bootstrap_model = self._bootstrap_from_fused(compute_model)
        return CoverageModel(
            xr.Dataset(
                dict(raw_model=raw_model, bootstrap_model=bootstrap_model),
                dict(factor=factors.factor_names),
            )
        )

    def predict(
        self,
        psicov: PsiCoverage,
        factors: Factors,
        base_psicov: Optional[PsiCoverage] = None,
        **update_attrs,
    ) -> "PsiCoverage":
        """Create updated :class:`PsiCoverage` corresponding to model and factors

        Parameters
        ----------
        psicov: PsiCoverage
        factors: Factors
        base_psicov: Optional[PsiCoverage]
            If specified, use this PsiCoverage object instead of psicov to add
            updated values to.
            Generally, we want to update `psicov`. But we don't want Dask to
            think it needs to rechunk more than necessary.
            It is most efficient to make predictions with:
            - `psicov` loaded with ec_idx_nchunks=None
            - `base_psicov` loaded with ec_idx_nchunks=1.
        **update_attrs
            Additional attributes to add to resulting PsiCoverage

        """
        # chunks to use
        prefix_chunks = psicov.raw_psi.chunksizes["prefix"]
        ec_idx_chunks = psicov.raw_psi.chunksizes["ec_idx"]
        # offsets on ec_idx
        offsets = psicov.lsv_offsets.astype(np.int64)
        # select factors
        try:
            factors = factors[psicov.prefixes]
        except KeyError as err:
            missing = set(psicov.prefixes) - set(factors.prefixes)
            raise ValueError(f"Missing input factors for prefixes: {missing}") from err
        factors_arr = factors[psicov.prefixes].factors.chunk({"prefix": prefix_chunks})
        # update raw/bootstrap separately
        updated = dict()
        model_pairs = {
            "raw_psi": (
                psicov.raw_psi,
                self.raw_model.chunk({"ec_idx": ec_idx_chunks}),
            ),
            "bootstrap_psi": (
                psicov.bootstrap_psi,
                self.bootstrap_model.chunk({"ec_idx": ec_idx_chunks}),
            ),
        }
        for key, (original_psi, model) in model_pairs.items():
            updated[key] = (
                original_psi.pipe(
                    mc.correct_with_model,
                    model,
                    factors_arr,
                    dim_prefix="prefix",
                    dim_factor="factor",
                )
                .pipe(
                    lambda unnormalized_psi: xr.apply_ufunc(
                        clip_and_normalize_strict,
                        unnormalized_psi.chunk({"ec_idx": -1}),
                        offsets,
                        input_core_dims=[["ec_idx"], ["e_offsets_idx"]],
                        output_core_dims=[["ec_idx"]],
                        dask="allowed",
                    )
                )
                .chunk({"ec_idx": ec_idx_chunks})
                .pipe(
                    lambda x: x.where(
                        x.notnull() & psicov.event_passed,
                        original_psi,  # noqa: B023
                    )
                )
            )
        # offsets for events for normalization
        if base_psicov is None:
            base_psicov = psicov
        return base_psicov.updated(
            bootstrap_psi=updated["bootstrap_psi"],
            raw_psi=updated["raw_psi"],
            **update_attrs,
        )
