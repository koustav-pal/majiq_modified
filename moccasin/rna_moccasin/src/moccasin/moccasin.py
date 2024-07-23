"""
moccasin.py

MOCCASIN core library

Author: Joseph K Aicher
"""

import itertools
import logging
import warnings
from pathlib import Path
from typing import Collection, Final, Optional, Union

import dask.array as da
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
from dask.distributed import default_client, progress

import moccasin.constants as constants
from moccasin.internals import infer_params as _infer_params
from moccasin.logger import get_logger


def persist_with_progress(x):
    try:
        default_client()  # raises ValueError if not using distributed cluster
    except ValueError:
        with ProgressBar():
            return x.persist()
    # we have a distributed client
    x = x.persist()  # start computation asynchronously
    # set up progress bar on running computation
    if isinstance(x, xr.Dataset):
        progress(*(v.data for v in x.variables.values() if v.chunks))
    elif isinstance(x, xr.DataArray) and x.chunks:
        progress(x.data)
    else:
        # assume it's a Dask object
        progress(x)
    return x


def clip_and_renormalize_psi(
    psi: xr.DataArray,
    offsets: Union[np.ndarray, xr.DataArray],
    dim_ecidx: str = constants.DEFAULT_DIM_ECIDX,
) -> xr.DataArray:
    """Clip values of PSI, renormalize each event to add up to 1 (nan if all 0)

    Parameters
    ----------
    psi: xr.DataArray
        Values to be clipped/normalized over required named dimension
        `dim_ecidx` relative to provided offsets
    offsets: Union[np.ndarray, xr.DataArray]
        Start/end offsets of slices to be normalized. offsets[0] = 0,
        offsets[-1] should equal length of `dim_ecidx` in psi.
        Normalization happens with respect to sum(psi[offsets[i]:offsets[i+1]])
    dim_ecidx: str
        Name of dimension in psi that normalization should happen against
    """
    if dim_ecidx not in psi.dims:
        raise ValueError(f"psi does not have required dim {dim_ecidx} ({psi.dims = })")
    # get offsets as ndarray
    offsets_arr = np.array(offsets, copy=False)
    # get function for sum over values between offsets over core dimension
    total_coverage_fn = _get_total_coverage_function(offsets_arr)
    # clip psi to 0, make ec_idx not a core dimension
    psi = psi.clip(min=0).chunk({dim_ecidx: -1})
    # normalization constant for psi
    psi_Z = xr.apply_ufunc(
        total_coverage_fn,
        psi,
        input_core_dims=[[dim_ecidx]],
        output_core_dims=[[dim_ecidx]],
        output_dtypes=(psi.dtype,),
        dask="parallelized",
    )
    # renormalize (nan mask if zero)
    return psi / psi_Z.where(psi_Z > 0)


def infer_model_params(
    uncorrected: xr.DataArray,
    passed: xr.DataArray,
    factors: xr.DataArray,
    complete: bool = False,
    dim_prefix: str = constants.DEFAULT_DIM_PREFIX,
    dim_factor: str = constants.DEFAULT_DIM_FACTOR,
) -> xr.DataArray:
    """get coefficients for OLS on uncorrected data given factors

    Parameters
    ----------
    uncorrected: xr.DataArray
        Observed values per prefix. NaN values are propagated regardless of
        passed vs not.
        Dimensions: (`dim_prefix`, ...)
    passed: xr.DataArray
        Boolean indicator if observation passed and used for modeling or
        ignored in modeling.
    factors: xr.DataArray
        Confounding and non-confounding factors which are used to predict the
        values in uncorrected. Must have all prefixes found in uncorrected
        Dimensions: (`dim_prefix`, `dim_factor`)
    complete: bool
        Indicates if we should handle missing data in uncorrected. If the data
        are complete, can be much faster doing this
    dim_prefix: str
        Name of dimension corresponding to experiments in uncorrected, factors
    dim_factor: str
        Name of dimension corresponding to factors in factors

    Returns
    -------
    xr.DataArray
        Coefficients such that xr.dot(factors, coeff, dim=dim_factor) is the OLS
        estimator of uncorrected given observed data
        Dimensions: (..., `dim_factor`)
    """
    # check that we have required dimensions
    if dim_prefix not in uncorrected.dims:
        raise ValueError(
            f"uncorrected does not have required dim {dim_prefix} ({uncorrected.dims = })"
        )
    if dim_prefix not in factors.dims:
        raise ValueError(
            f"factors does not have required dim {dim_prefix} ({factors.dims = })"
        )
    if dim_factor not in factors.dims:
        raise ValueError(
            f"factors does not have required dim {dim_factor} ({factors.dims = })"
        )
    DIM_EXTRA = "_extra_core_dims"  # put all dimensions together if have extra
    if DIM_EXTRA in (dim_factor, dim_prefix):
        raise ValueError(
            f"dim_factor, dim_prefix cannot be set to reserved name {DIM_EXTRA}"
        )
    # get factors only for prefixes in uncorrected, make same dtype as uncorrected
    factors = factors.sel(prefix=uncorrected[dim_prefix]).astype(uncorrected.dtype)
    # get original names with indexes from uncorrected
    uncorrected_indexes = set(uncorrected.indexes)
    extra_core_dims = list(
        set(uncorrected.dims) - (set(factors.dims) | set(passed.dims))
    )
    if extra_core_dims:
        uncorrected = uncorrected.stack({DIM_EXTRA: extra_core_dims}).chunk(
            {DIM_EXTRA: -1}
        )
    else:
        uncorrected = uncorrected.expand_dims(dim=DIM_EXTRA)
    uncorrected = uncorrected.chunk({dim_prefix: -1})
    passed = passed.chunk({dim_prefix: -1})
    params = xr.apply_ufunc(
        _infer_params,
        uncorrected,
        factors,
        passed,
        input_core_dims=[
            [DIM_EXTRA, dim_prefix],
            [dim_factor, dim_prefix],
            [dim_prefix],
        ],
        output_core_dims=[[DIM_EXTRA, dim_factor]],
        dask="allowed",
        output_dtypes=[uncorrected.dtype],
    )
    if extra_core_dims:
        # unstack: get parameters back in original dimensions
        params = params.unstack(DIM_EXTRA)
    else:
        params = params.squeeze(dim=DIM_EXTRA, drop=True)
    # don't keep dimension labels that didn't already exist
    params = params.drop_vars(
        [
            x
            for x in params.coords
            if not (x in uncorrected_indexes or x in factors.indexes)
        ]
    )
    return params


def correct_with_model(
    uncorrected: xr.DataArray,
    params: xr.DataArray,
    factors: xr.DataArray,
    dim_prefix: str = constants.DEFAULT_DIM_PREFIX,
    dim_factor: str = constants.DEFAULT_DIM_FACTOR,
) -> xr.DataArray:
    """add residuals from confounded model to unconfounded predictions

    Notes
    -----
    These are corrected values under linear model. Does not include any
    transformations (e.g. clipping, log-transformation, renormalization)

    Parameters
    ----------
    uncorrected: xr.DataArray
        Observed values per prefix. NaN values indicate missing data.
        Dimensions: (?`dim_prefix`, ...)
    params: xr.DataArray
        Model parameters from infer_model_params()
        Dimensions: (..., `dim_factor`)
    factors: xr.DataArray
        Confounding and non-confounding factors which are used to predict the
        values in uncorrected. Must have all prefixes found in uncorrected if
        uncorrected has prefixes (if uncorrected does not have prefix, then
        factors must not either).
        Coordinate "confounding" is boolean mask over factors indicating if
        factor is confounding vs unconfounding
        Dimensions: (?`dim_prefix`, `dim_factor`)
        Coordinates: confounding (dim: `dim_factor`)
    dim_prefix: str
        Name of dimension corresponding to experiments
    dim_factor: str
        Name of dimension corresponding to factors

    Returns
    -------
    xr.DataArray
        Corrected values. Residuals from predictions with all factors (vs true
        values) added to predictions with nonconfounding factors only
    """
    # check for required dimensions
    if dim_factor not in factors.dims:
        raise ValueError(
            f"factors does not have required dim {dim_factor} ({factors.dims = })"
        )
    if dim_factor not in params.dims:
        raise ValueError(
            f"params does not have required dim {dim_factor} ({params.dims = })"
        )
    # make sure prefixes are aligned appropriately
    if dim_prefix in uncorrected.dims:
        # if uncorrected has prefixes, we expect it in factors, require match
        factors = factors.sel({dim_prefix: uncorrected[dim_prefix]})
    # if they are all nonconfounding, correction should yield uncorrected:
    if not factors["confounding"].any():
        return uncorrected

    # make sure we don't make uncorrected higher precision on accident
    factors = factors.astype(uncorrected.dtype)  # same type as uncorrected

    # theoretically, the model can be written in pseudocode as:
    # predictions = xr.dot(factors, params, dims=dim_factor)  # model with all factors
    # residuals = uncorrected - predictions  # residuals from prediction
    # unconfounded_predictions = xr.dot(factors.where(not confounding, 0), params, dims=dim_factor)
    # return unconfounded_predictions + residuals

    # that is: get residuals from full model, then add it back to predictions
    # from model without confounding values.
    # in the linear model, by associative/distributive properties, this can be
    # simplified into one matrix dot product:
    confounding_idx = np.where(factors["confounding"].values)[0]
    return uncorrected - xr.dot(
        params.isel({dim_factor: confounding_idx}),
        factors.isel({dim_factor: confounding_idx}),
        dim=dim_factor,
    )


class ModelUnknownConfounders(object):
    """Model for how to use observations and factors to find unknown confounders"""

    def __init__(
        self,
        original_ecidx: xr.DataArray,
        model_params: xr.DataArray,
        singular_values: xr.DataArray,
        ec_vectors: xr.DataArray,
        total_variance: float,
        dim_factor: str,
        dim_ecidx: str,
        default_num_factors: Optional[int] = None,
    ):
        """Model of unknown confounders with specified parameters

        Parameters
        ----------
        original_ecidx: xr.DataArray
            indexes into ec_idx from original event connections
            (dimension: top_`dim_ecidx`)
        model_params: xr.DataArray
            coefficients for linear model to get residuals on which unknown
            confounders are computed
            (dimension: top_`dim_ecidx`, `dim_factor`)
        singular_values: xr.DataArray
            singular values for unknown confounders found from training data
            (dimension: new_`dim_factor`)
        ec_vectors: xr.DataArray
            right-singular vectors over top_ec_idx used to project to new
            factors (scaled appropriately by singular values)
            (dimension: new_`dim_factor`, top_`dim_ecidx`)
        total_variance: float
            sum of variances scaled by the number of experiments used to train.
            This allows us to calculate the variance explained by the singular
            values
        dim_factor: str
            Base name of dimension corresponding to factors
        dim_ecidx: str
            Base name of dimension corresponding to event connections
        default_num_factors: Optional[int]
            Default number of factors to use when inferring new factors from
            the model. When None, use all possible new factors in the model.
        """
        # model values
        self.original_ecidx: Final[xr.DataArray] = original_ecidx
        self.model_params: Final[xr.DataArray] = model_params
        self.singular_values: Final[xr.DataArray] = singular_values
        self.ec_vectors: Final[xr.DataArray] = ec_vectors
        # total variance, used to calculate variance explained
        self.total_variance: Final[float] = total_variance
        # base names for factor/ecidx dimensions
        self.dim_factor: Final[str] = dim_factor
        self.dim_ecidx: Final[str] = dim_ecidx
        # TODO: check that dimensions are as expected
        dim_top_ecidx = self.get_dim_top_ecidx(self.dim_ecidx)
        dim_new_factor = self.get_dim_new_factor(self.dim_factor)
        if dim_top_ecidx not in original_ecidx.dims:
            raise ValueError(
                f"original_ecidx does not have required dim {dim_top_ecidx}"
                f" ({original_ecidx.dims = })"
            )
        if dim_top_ecidx not in model_params.dims:
            raise ValueError(
                f"model_params does not have required dim {dim_top_ecidx}"
                f" ({model_params.dims = })"
            )
        if dim_factor not in model_params.dims:
            raise ValueError(
                f"model_params does not have required dim {dim_factor}"
                f" ({model_params.dims = })"
            )
        if dim_new_factor not in singular_values.dims:
            raise ValueError(
                f"singular_values does not have required dim {dim_new_factor}"
                f" ({singular_values.dims = })"
            )
        if dim_new_factor not in ec_vectors.dims:
            raise ValueError(
                f"ec_vectors does not have required dim {dim_new_factor}"
                f" ({ec_vectors.dims = })"
            )
        if dim_top_ecidx not in ec_vectors.dims:
            raise ValueError(
                f"ec_vectors does not have required dim {dim_top_ecidx}"
                f" ({ec_vectors.dims = })"
            )
        # default num new factors must be valid
        self.default_num_factors: Final[int] = (
            default_num_factors or self.max_num_factors
        )
        if self.default_num_factors > self.max_num_factors:
            raise ValueError(
                f"{default_num_factors = } is greater than available number of"
                f" singular values ({self.max_num_factors})"
            )
        elif self.default_num_factors < 0:
            raise ValueError(f"{default_num_factors = } must be non-negative")
        return

    @property
    def max_num_factors(self) -> int:
        """Maximum number of new factors to infer limited by stored components"""
        return self.singular_values.sizes[self.dim_new_factor]

    @staticmethod
    def get_dim_top_ecidx(dim_ecidx: str) -> str:
        """translate dim_ecidx to dimension for top ecidx"""
        return f"top_{dim_ecidx}"

    @staticmethod
    def get_dim_new_factor(dim_factor: str) -> str:
        """translate dim_factor to dimension for new factors"""
        return f"new_{dim_factor}"

    @property
    def dim_top_ecidx(self) -> str:
        """translate dim_ecidx to dimension for top ecidx"""
        return self.get_dim_top_ecidx(self.dim_ecidx)

    @property
    def dim_new_factor(self) -> str:
        """translate dim_factor to dimension for new factors"""
        return self.get_dim_new_factor(self.dim_factor)

    def to_zarr(self, output: Union[str, Path]) -> None:
        """Save model parameters fo file"""
        xr.Dataset(
            {
                "original_ecidx": self.original_ecidx,
                "model_params": self.model_params,
                "singular_values": self.singular_values,
                "ec_vectors": self.ec_vectors,
            },
            {},
            {
                "total_variance": self.total_variance,
                "dim_factor": self.dim_factor,
                "dim_ecidx": self.dim_ecidx,
                "default_num_factors": self.default_num_factors,
            },
        ).to_zarr(output, mode="w", consolidated=True)
        return

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "ModelUnknownConfounders":
        """Load model parameters from file"""
        with xr.open_zarr(path) as df:
            df.load()
            return ModelUnknownConfounders(
                original_ecidx=df.original_ecidx,
                model_params=df.model_params,
                singular_values=df.singular_values,
                ec_vectors=df.ec_vectors,
                total_variance=df.attrs["total_variance"],
                dim_factor=df.attrs.get("dim_factor", constants.DEFAULT_DIM_FACTOR),
                dim_ecidx=df.attrs.get("dim_ecidx", constants.DEFAULT_DIM_ECIDX),
                default_num_factors=df.attrs.get("default_num_factors", None),
            )

    @classmethod
    def train(
        cls,
        uncorrected: xr.DataArray,
        passed: xr.DataArray,
        offsets: np.ndarray,
        factors: xr.DataArray,
        default_num_factors: int = constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
        max_new_factors: Optional[int] = None,
        max_events: int = constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        dim_prefix: str = constants.DEFAULT_DIM_PREFIX,
        dim_factor: str = constants.DEFAULT_DIM_FACTOR,
        dim_ecidx: str = constants.DEFAULT_DIM_ECIDX,
        log: Optional[logging.Logger] = None,
        show_progress: bool = False,
    ) -> "ModelUnknownConfounders":
        """Learn model for unknown confounders using residuals from OLS model

        Parameters
        ----------
        uncorrected: xr.DataArray (dims: `dim_prefix`, ..., `dim_ecidx`)
            data from which variations are used to identify confounding factors
        passed: xr.DataArray (dims: `dim_prefix`, ..., `dim_ecidx`)
            indicator as to whether values in uncorrected passed and should be
            included in modeling
        offsets: np.ndarray
            offsets into ec_idx describing different events (only use 0 or 1
            connection per event in model)
        factors: xr.DataArray (dims: `dim_prefix`, `dim_factor`)
            confounding/nonconfounding factors used to predict uncorrected.
            Residuals from OLS model used for SVD.
            Note that "confounding" coordinate, if present, is ignored.
        default_num_factors: int
            default number of factors that the model will output when not
            otherwise specified in the `predict` method.
        max_new_factors: Optional[int]
            Deprecated, the model will always solve for parameters covering the
            maximum possible number of new factors.
            Replaced by `default_num_factors`, which determines default number
            of factors to output (when not otherwise specified in the `predict`
            method).
            If set, triggers deprecation warning and overrides
            `default_num_factors` for backwards compatibility.
        max_events: int
            maximum number of events used to create the model
        dim_prefix: str
            Name of dimension corresponding to experiments
        dim_factor: str
            Name of dimension corresponding to factors
        dim_ecidx: str
            Name of dimension corresponding to event connections
        show_progress: bool
            Show progress bar for computations, assuming dask distributed
            scheduler

        Returns
        -------
        ModelUnknownConfounders
            Model picks top events that are fully quantified and have highest
            variance.
            Model uses factors to build model predicting uncorrected.
            Residuals are used as input into SVD, and appropriate parameters
            saved to be able to produce new unknown confounders for any future
            sample
        """
        if log is None:
            log = get_logger()
        if max_new_factors is not None:
            warnings.warn(
                "max_new_factors is deprecated, use default_num_factors instead",
                DeprecationWarning,
                stacklevel=2,
            )
            default_num_factors = max_new_factors
        # check that inputs have required dimensions
        if dim_prefix not in uncorrected.dims:
            raise ValueError(
                f"uncorrected does not have required dim {dim_prefix}"
                f" ({uncorrected.dims = })"
            )
        if dim_ecidx not in uncorrected.dims:
            raise ValueError(
                f"uncorrected does not have required dim {dim_ecidx}"
                f" ({uncorrected.dims = })"
            )
        if dim_prefix not in passed.dims:
            raise ValueError(
                f"passed does not have required dim {dim_prefix}" f" ({passed.dims = })"
            )
        if dim_ecidx not in passed.dims:
            raise ValueError(
                f"passed does not have required dim {dim_ecidx}" f" ({passed.dims = })"
            )
        if dim_prefix not in factors.dims:
            raise ValueError(
                f"factors does not have required dim {dim_prefix}"
                f" ({factors.dims = })"
            )
        if dim_factor not in factors.dims:
            raise ValueError(
                f"factors does not have required dim {dim_factor}"
                f" ({factors.dims = })"
            )
        # get names of new dimensions used (top_ec_idx, new_factor)
        dim_top_ecidx = cls.get_dim_top_ecidx(dim_ecidx)
        dim_new_factor = cls.get_dim_new_factor(dim_factor)
        # load factors, make sure it's full rank
        factors = (
            factors.sel({dim_prefix: uncorrected[dim_prefix]})
            .astype(uncorrected.dtype)
            .load()
        )
        if (factors_rank := np.linalg.matrix_rank(factors.values)) < factors.sizes[
            dim_factor
        ]:
            raise ValueError(
                f"factors is not full rank ({factors_rank} < {factors.sizes[dim_factor]}"
            )
        # get indexes for connections for top LSVs
        uncorrected_medians = cls._median_extra_dims(
            uncorrected.where(passed), keep_dims=(dim_prefix, dim_ecidx)
        )
        log.info(
            "We expect likely NaN slice, invalid value in true_divide warnings"
            " in this next step (https://github.com/dask/dask/issues/3245)."
            " They can be safely ignored"
        )
        # get variance across samples of psi (median over bootstraps, etc.)
        original_ecidx_ds = xr.Dataset(
            {
                "variance": uncorrected_medians.var(dim_prefix),
                "complete": uncorrected_medians.count(dim_prefix)
                == uncorrected.sizes[dim_prefix],
            },
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "invalid value encountered in .*divide", RuntimeWarning
            )
            if show_progress:
                original_ecidx_ds = persist_with_progress(original_ecidx_ds)
            original_ecidx_ds = original_ecidx_ds.load()
        log.info("We do not expect warnings after this point")
        original_ecidx = (
            original_ecidx_ds
            # compute into in-memory pandas dataframe
            .to_dataframe()
            # get indexes with maximum variance, only 1 per LSV
            .assign(lsv_idx=np.repeat(np.arange(len(offsets) - 1), np.diff(offsets)))
            # keep only complete indexes, get top variances, 1 per LSV
            .loc[lambda df: df["complete"]]
            .sort_values("variance", ascending=False)
            .drop_duplicates("lsv_idx", keep="first")
            .head(max_events)
            # make appropriate DataArray of indexes for this
            .pipe(
                lambda x: xr.DataArray(
                    x.index.sort_values().values, dims=[dim_top_ecidx]
                )
            )
        )
        # get data with which we will train the model
        # TODO maybe want to persist this in memory if it isn't too big
        x = cls._median_extra_dims(
            uncorrected.isel({dim_ecidx: original_ecidx}),
            keep_dims=(dim_prefix, dim_top_ecidx),
        )
        # by definition of x, we know that all are passed
        x_passed = xr.DataArray(True).expand_dims(
            {dim_prefix: x[dim_prefix], dim_top_ecidx: x[dim_top_ecidx]}  # type: ignore[arg-type]
        )
        model_params = infer_model_params(
            x,
            x_passed,
            factors,
            complete=True,
            dim_prefix=dim_prefix,
            dim_factor=dim_factor,
        )
        x_residuals = x - xr.dot(factors, model_params, dim=dim_factor)
        # perform sparse SVD on x_residuals
        prefix_vectors, singular_values, ec_vectors = da.linalg.svd(
            x_residuals.chunk({dim_top_ecidx: -1})
            .transpose(dim_prefix, dim_top_ecidx)
            .data
        )
        # load model parameters into memory
        local = xr.Dataset(
            {
                "model_params": model_params,
                "prefix_vectors": xr.DataArray(
                    prefix_vectors,
                    dims=[dim_prefix, dim_new_factor],
                ),
                "singular_values": xr.DataArray(singular_values, dims=[dim_new_factor]),
                "ec_vectors": xr.DataArray(
                    ec_vectors,
                    dims=[dim_new_factor, dim_top_ecidx],
                ),
                "total_variance": xr.dot(x_residuals, x_residuals),
            },
            {
                dim_prefix: factors[dim_prefix],
                dim_new_factor: list(
                    itertools.islice(
                        # make sure there are no conflicts in names of unknown factors
                        (
                            name
                            for name in (f"_U_{n}" for n in itertools.count(1))
                            if name not in factors[dim_factor]
                        ),
                        len(singular_values),
                    )
                ),
                "confounding": (
                    dim_new_factor,
                    [True for _ in range(len(singular_values))],
                ),
            },
        )
        if show_progress:
            local = persist_with_progress(local)
        local = local.load()
        # only keep as many new factors that keep model matrix having full
        # column rank
        combined_factors = np.empty(
            (
                factors.sizes[dim_prefix],
                factors.sizes[dim_factor] + len(singular_values),
            ),
            dtype=factors.dtype,
        )
        combined_factors[:, : factors.sizes[dim_factor]] = factors.transpose(
            dim_prefix, dim_factor
        ).values
        combined_factors[
            :, factors.sizes[dim_factor] :
        ] = local.prefix_vectors.transpose(dim_prefix, dim_new_factor).values
        keep_n: int  # we know we can keep at least this many new factors (starting at 0)
        for keep_n in range(len(singular_values)):
            # try adding 1 more after the keep_n we know works, still full rank?
            total = factors.sizes[dim_factor] + keep_n + 1
            if np.linalg.matrix_rank(combined_factors[:, :total]) < total:
                # combined factors now singular so keep_n is maximum new factors
                break
        else:
            # all of combined_factors was still full rank
            keep_n = len(singular_values)

        return ModelUnknownConfounders(
            original_ecidx=original_ecidx,
            model_params=local.model_params,
            singular_values=local.singular_values.isel(new_factor=slice(keep_n)),
            ec_vectors=local.ec_vectors.isel(new_factor=slice(keep_n)),
            total_variance=local.total_variance.values[()],
            dim_factor=dim_factor,
            dim_ecidx=dim_ecidx,
            default_num_factors=min(default_num_factors, keep_n),
        )

    @property
    def explained_variance(self) -> xr.DataArray:
        """proportion of training data residuals variance explained by each new factor"""
        return (
            (self.singular_values * self.singular_values) / self.total_variance
        ).rename("explained_variance")

    @staticmethod
    def _median_extra_dims(x: xr.DataArray, keep_dims: Collection[str]) -> xr.DataArray:
        """collapse any dimensions that aren't prefix or ec_idx by median"""
        extra_dims = [x for x in x.dims if x not in keep_dims]
        return x.median(extra_dims)

    def predict(
        self,
        uncorrected: xr.DataArray,
        passed: xr.DataArray,
        factors: xr.DataArray,
        dim_prefix: str = constants.DEFAULT_DIM_PREFIX,
        num_factors: Optional[int] = None,
    ) -> xr.DataArray:
        """get unknown confounders using observed data/factors

        Parameters
        ----------
        uncorrected: xr.DataArray
            Dimensions: (?dim_prefix, ..., self.dim_ecidx)
        passed: xr.DataArray
            Dimensions: (?dim_prefix, ..., self.dim_ecidx)
        factors: xr.DataArray
            Dimensions: (?dim_prefix, ..., self.dim_factor)
        dim_prefix: str
            dimension where prefix information (if present) is found
        num_factors: Optional[int]
            If specified, number of unknown confounding factors to infer.
            Must be >= 0 and <= `self.max_num_factors`.
            If None, use `self.default_num_factors`

        Notes
        -----
        Unpassed values have their residuals imputed to 0 (by
        definition experiments used in training must be present, but not true
        for held out data)
        """
        # check for required dimensions
        if self.dim_ecidx not in uncorrected.dims:
            raise ValueError(
                f"uncorrected does not have required dim {self.dim_ecidx}"
                f" ({uncorrected.dims = })"
            )
        if self.dim_ecidx not in passed.dims:
            raise ValueError(
                f"passed does not have required dim {self.dim_ecidx}"
                f" ({passed.dims = })"
            )
        if self.dim_factor not in factors.dims:
            raise ValueError(
                f"factors does not have required dim{self.dim_factor}"
                f" ({factors.dims = })"
            )
        # make sure prefixes aligned appropriately
        if dim_prefix in uncorrected.dims:
            # if uncorrected has prefixes, we expect it in factors, require match
            factors = factors.sel({dim_prefix: uncorrected[dim_prefix]})
        # make sure factors doesn't cause us to make uncorrected higher precision
        factors = factors.astype(uncorrected.dtype)
        # get subset of uncorrected to use
        x = uncorrected.where(passed)  # mask unpassed values
        x = x.isel({self.dim_ecidx: self.original_ecidx})  # select subset
        x = self._median_extra_dims(x, keep_dims=(dim_prefix, self.dim_top_ecidx))
        # get residuals, imputing 0 for missing values possible in held-out data
        residuals = x - xr.dot(factors, self.model_params, dim=self.dim_factor)
        residuals = residuals.where(residuals.notnull(), 0)
        # get singular values/ec_vectors to project with
        if num_factors is None:
            num_factors = self.default_num_factors
        if num_factors > self.max_num_factors:
            raise ValueError(
                f"{num_factors = } is greater than available number of"
                f" singular values ({self.max_num_factors})"
            )
        elif num_factors < 0:
            raise ValueError(f"{num_factors = } must be non-negative")
        num_factors_slice = {self.dim_new_factor: slice(num_factors)}
        # project residuals to new factors
        return xr.dot(
            residuals,
            np.reciprocal(self.singular_values.isel(num_factors_slice)),
            self.ec_vectors.isel(num_factors_slice),
            dim=self.dim_top_ecidx,
        )


def _get_total_coverage_function(offsets: np.ndarray):
    """get function that summarizes totals defined by offsets"""
    offsets = offsets.astype(np.int64)
    event_size = np.diff(offsets)

    def result(x: np.ndarray) -> np.ndarray:
        """get total coverage within each event defined by offsets"""
        if x.shape[-1] != offsets[-1]:
            raise ValueError(
                "input array core dim size does not match offsets"
                f" ({x.shape = }, {offsets[-1] = })"
            )
        return np.repeat(np.add.reduceat(x, offsets[:-1], axis=-1), event_size, axis=-1)

    return result
