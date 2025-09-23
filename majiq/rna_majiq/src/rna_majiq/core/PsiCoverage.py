"""
PsiCoverage.py

PSI and total coverage (raw and bootstrapped). Converted to/from EventsCoverage.
This allows simplification of the workflow with MOCCASIN bootstrap correction
and more readily parallelized analysis of arbitrarily many files by handling
dependences between junctions.

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Hashable,
    List,
    Optional,
    Sequence,
    Union,
    cast,
)

import dask.array as da
import numpy as np
import numpy.typing as npt
import xarray as xr
from dask.delayed import Delayed
from moccasin.moccasin import persist_with_progress

import rna_majiq._offsets as _offsets
import rna_majiq.constants as constants
from rna_majiq._version import version as nm_version
from rna_majiq.logger import get_logger

from .Events import Events
from .EventsCoverage import EventsCoverage
from .MixinPsiInference import MixinBootstrapPsi, MixinPsiOverPrefixes
from .SJExperiment import SJExperiment

if TYPE_CHECKING:
    from .PsiGroup import PsiGroup


class PsiCoverage(MixinBootstrapPsi, MixinPsiOverPrefixes):
    """Summarized raw and bootstrap coverage over LSVs for one or more experiments.

    Summarized raw and bootstrap coverage over LSVs for one or more
    experiments as input for quantification.
    Coverage is a total readrate over all bins, excluding stacks and after any
    preceding batch correction steps, ready for quantification.
    Per-experiment coverage stored independently over "prefix" dimension, where
    prefixes originate as the prefix from BAM file names
    (i.e. foo/experiment1.bam -> experiment1).
    Coverage is accompanied by boolean array indicating whether an event is
    "passed for quantification" for each experiment
    (:attr:`PsiCoverage.passed`).

    Provides functionality for combining and summarizing over experiments and
    multiple :class:`PsiCoverage` objects.
    Functions and attributes enable computation of PSI posterior statistics
    under MAJIQ models for splicing quantification.
    Computations are performed over xarray objects.
    When loading :class:`PsiCoverage` from Zarr files, data/computations
    will be loaded/performed lazily using Dask.
    Testing of these computations have been performed over local clusters using
    threads rather than processes (expensive computations generally release the
    GIL).

    Generally, for point estimates of location, quantification with raw coverage
    should be preferred, as bootstrap estimates converge very closely to raw
    estimates as the number of bootstrap replicates becomes large.
    For estimates of variability, quantification with bootstrap coverage should
    be used to account for additional per-bin readrate variability that isn't
    fully captured by the Bayesian model on its own.

    Underlying coverage is stored as the total number of reads over the event
    and the proportion of reads per intron/junction.
    This requires twice the uncompressed memory vs the number of reads per
    intron/junction, but permits easier lazy computation with Dask over large
    datasets.

    Parameters
    ----------
    df: xr.Dataset
        Required variables/coordinates as in `EXPECTED_VARIABLES`
    events: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events

    See Also
    --------
    PsiCoverage.from_sj_lsvs : Create :class:`PsiCoverage` from :class:`SJExperiment` and :class:`Events`
    PsiCoverage.from_events_coverage : Create :class:`PsiCoverage` from :class:`EventsCoverage`
    PsiCoverage.from_zarr : Load :class:`PsiCoverage` from one or more Zarr files
    PsiCoverage.updated : Create updated :class:`PsiCoverage` with updated arrays
    PsiCoverage.sum : Summed :class:`PsiCoverage` over current prefixes
    PsiCoverage.mask_events : Create updated :class:`PsiCoverage` passing only specified events
    PsiCoverage.__getitem__ : Get :class:`PsiCoverage` for subset of prefixes
    """

    EXPECTED_VARIABLES = {
        "event_passed": ("ec_idx", "prefix"),
        "raw_total": ("ec_idx", "prefix"),
        "raw_psi": ("ec_idx", "prefix"),
        "bootstrap_total": ("ec_idx", "prefix", "bootstrap_replicate"),
        "bootstrap_psi": ("ec_idx", "prefix", "bootstrap_replicate"),
        "prefix_total": ("prefix",),
        "prefix": ("prefix",),
    }

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Initialize :class:`PsiCoverage` with specified xarray datasets

        Parameters
        ----------
        df: xr.Dataset
            Required variables/coordinates as in `EXPECTED_VARIABLES`
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        """
        # save PSI information and events
        MixinPsiOverPrefixes.__init__(self, df, events)
        return

    def group(
        self,
        save_zarr: Optional[Union[str, Path]] = None,
        tmp_base_dir: Optional[Union[str, Path]] = None,
        ec_idx_nchunks: Optional[int] = 1,
        show_progress: bool = False,
    ) -> "PsiGroup":
        """Create :class:`PsiGroup` from coverage in `self`

        Parameters
        ----------
        save_zarr: Optional[Union[str, Path]]
            If specified, save rechunked data permanently to specified path.
            Otherwise (default), use a temporary directory.
        tmp_base_dir: Optional[Union[str, Path]]
            Set base directory for temporary directories with rechunked data.
            Use if default from `tempfile` is inappropriate.
        ec_idx_nchunks: Optional[int]
            How many psigroup chunks to load for ec_idx at a time when
            computing with dask
        show_progress: bool
            Show progressbar when rechunking data.

        Returns
        -------
        PsiGroup
            Rechunked raw and approximate posterior distribution of PSI
            parameters to prioritize computations over prefixes
        """
        from .PsiGroup import PsiGroup

        return PsiGroup.from_psicov(
            self,
            save_zarr=save_zarr,
            tmp_base_dir=tmp_base_dir,
            ec_idx_nchunks=ec_idx_nchunks,
            show_progress=show_progress,
        )

    @property
    def num_bootstraps(self) -> int:
        """Number of bootstrap replicates used for bootstraped coverage estimates"""
        return self.df.sizes["bootstrap_replicate"]

    @cached_property
    def event_passed(self) -> xr.DataArray:
        """array(prefix, ec_idx) indicating if event passed"""
        return self.df["event_passed"].astype(bool)

    @property
    def prefix_total(self) -> xr.DataArray:
        """array(prefix) of total number of reads over entire experiment"""
        return cast(xr.DataArray, self.df["prefix_total"].reset_coords(drop=True))

    @property
    def raw_total(self) -> xr.DataArray:
        """array(prefix, ec_idx) raw total reads over event

        py:class:`xr.DataArray` (prefix, ec_idx) raw total reads over event
        (i.e. sum over all event connections per event)
        """
        return cast(xr.DataArray, self.df["raw_total"].reset_coords(drop=True))

    @property
    def raw_psi(self) -> xr.DataArray:
        """array(prefix, ec_idx) percentage of raw_total for connection

        py:class:`xr.DataArray` (prefix, ec_idx) percentage of raw_total for
        connection (maximum likelihood estimate of PSI over raw reads)
        """
        return cast(xr.DataArray, self.df["raw_psi"].reset_coords(drop=True))

    @property
    def bootstrap_total(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_total

        py:class:`xr.DataArray` (prefix, ec_idx, bootstrap_replicate)
        bootstrapped total reads over event (i.e. sum over all event
        connections per event)
        """
        return cast(xr.DataArray, self.df["bootstrap_total"].reset_coords(drop=True))

    @property
    def bootstrap_psi(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_psi

        py:class:`xr.DataArray` (prefix, ec_idx, bootstrap_replicate)
        percentage of bootstrap_total for connection (maximum likelihood
        estimate of PSI over raw reads)
        """
        return cast(xr.DataArray, self.df["bootstrap_psi"].reset_coords(drop=True))

    @cached_property
    def raw_coverage(self) -> xr.DataArray:
        """array(prefix, ec_idx) coverage for individual connection (psi * total)"""
        return self.raw_psi * self.raw_total

    @cached_property
    def bootstrap_coverage(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_coverage"""
        return self.bootstrap_psi * self.bootstrap_total

    @cached_property
    def alpha_prior(self) -> xr.DataArray:
        """array(ec_idx) alpha parameter of prior distribution on PSI for connection"""
        return 1 / self.event_size.astype(self.raw_psi.dtype)

    @cached_property
    def beta_prior(self) -> xr.DataArray:
        """array(ec_idx) beta parameter of prior distribution on PSI for connection"""
        return 1 - self.alpha_prior

    @cached_property
    def raw_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx) alpha parameter of raw posterior"""
        return (self.raw_coverage + self.alpha_prior).where(self.event_passed)

    @cached_property
    def bootstrap_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) alpha parameter of bootstrapped posterior"""
        return (self.bootstrap_coverage + self.alpha_prior).where(self.event_passed)

    @cached_property
    def raw_alpha_plus_beta(self) -> xr.DataArray:
        return 1 + self.raw_total

    @cached_property
    def raw_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx) beta parameter of raw posterior"""
        return self.raw_alpha_plus_beta - self.raw_alpha

    @cached_property
    def bootstrap_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) beta parameter of bootstrapped posterior"""
        return 1 + self.bootstrap_total - self.bootstrap_alpha

    @classmethod
    def from_events_coverage(
        cls,
        events_coverage: EventsCoverage,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        prefix_total: float = np.nan,
    ) -> "PsiCoverage":
        """Create :py:class:`PsiCoverage` from :py:class:`EventsCoverage`

        Parameters
        ----------
        events_coverage: EventsCoverage
            Read coverage per feature, raw and bootstrapped, along with raw
            number of bins with nonzero coverage
        minreads, minbins: float
            Quantifiability thresholds
        prefix_total: float
            If specified, the sum total number of reads over all junctions for
            the experiment

        Returns
        -------
        PsiCoverage
        """
        # get offsets as int (not uint)
        offsets: npt.NDArray[np.int64] = cast(
            npt.NDArray[np.int64],
            events_coverage.events._offsets.view(np.int64),
        )
        # get whether individual connection passes thresholds
        passed = (events_coverage.numreads >= minreads) & (
            events_coverage.numbins >= minbins
        )
        # get whether any connection in event passed, per connection
        event_passed = _offsets.offset_logical_or(passed, offsets)
        # get total coverage per event, per connection
        raw_total = _offsets.offsetsum(
            events_coverage.numreads, offsets, axes=[0, -1, 0]
        )
        bootstrap_total = _offsets.offsetsum(
            events_coverage.bootstraps, offsets, axes=[0, -1, 0]
        )
        # get psi per connection
        with np.errstate(divide="ignore", invalid="ignore"):
            raw_psi = np.where(raw_total > 0, events_coverage.numreads / raw_total, 0)
            bootstrap_psi = np.where(
                bootstrap_total > 0, events_coverage.bootstraps / bootstrap_total, 0
            )
        # return dataset with matched values
        return cls(
            xr.Dataset(
                data_vars=dict(
                    event_passed=("ec_idx", event_passed),
                    raw_total=("ec_idx", raw_total),
                    raw_psi=("ec_idx", raw_psi),
                    bootstrap_total=(
                        ("ec_idx", "bootstrap_replicate"),
                        bootstrap_total,
                    ),
                    bootstrap_psi=(("ec_idx", "bootstrap_replicate"), bootstrap_psi),
                    prefix_total=np.array(prefix_total, dtype=raw_total.dtype),
                ),
                attrs=dict(
                    minreads=minreads,
                    minbins=minbins,
                    bam_path=events_coverage.bam_path,
                    bam_version=events_coverage.bam_version,
                ),
            ).expand_dims(prefix=[events_coverage.prefix]),
            events_coverage.events.save_df,
        )

    @classmethod
    def from_sj_lsvs(
        cls,
        sj: SJExperiment,
        lsvs: Events,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> "PsiCoverage":
        """Create :class:`PsiCoverage` from :class:`SJExperiment` and :class:`Events`.

        Parameters
        ----------
        sj: SJExperiment
            Intron and junction coverage from an experiment
        lsvs: Events
            Events over which PsiCoverage will be defined
        minreads, minbins: float
            Quantifiability thresholds
        num_bootstraps: int
            The number of bootstrap replicates for bootstrapped estimates
        pvalue_threshold: float
            P-value threshold for removing stacks under leave-one-out Poisson
            model of per-bin read coverage, for both raw and bootstrapped
            coverage (Set to nonpositive value to skip stack detection)

        Returns
        -------
        PsiCoverage

        Notes
        -----
        The pvalue_threshold for stack removal is applied to both raw and
        bootstrapped coverage. This differs from the behavior in MAJIQ v2,
        where stack removal was only applied to bootstrapped coverage.
        In this sense "raw" coverage is only after stack detection.
        """
        lsv_coverage = EventsCoverage.from_events_and_sj(
            lsvs, sj, num_bootstraps=num_bootstraps, pvalue_threshold=pvalue_threshold
        )
        return PsiCoverage.from_events_coverage(
            lsv_coverage,
            minreads=minreads,
            minbins=minbins,
            prefix_total=sj.junctions.numreads(
                numstacks=sj.junctions.numstacks(pvalue_threshold=pvalue_threshold)
            ).sum(),
        )

    @classmethod
    def from_zarr(
            cls,
            path: Union[str, Path, List[Union[str, Path]]],
            ec_idx_nchunks: Optional[int] = None,
            prefix_nchunks: Optional[int] = 1,
            preload: Optional[bool] = False
    ) -> "PsiCoverage":
        """Load :py:class:`PsiCoverage` from one or more specified paths

        Load :py:class:`PsiCoverage` from one or more specified paths.
        Prefixes will be concatenated (overlapping prefixes will use values
        from the first file it is found in).

        Parameters
        ----------
        path: Union[str, Path, List[Union[str, Path]]]
            path or paths with PsiCoverage saved in zarr format
        ec_idx_nchunks, prefix_nchunks: Optional[int]
            Number of chunks on disk to group per chunk managed by dask.
            If None, load all chunks along the dimension together.
        preload: Optional[bool]
            If set to true, on-disk zarr array will be read into memory
            by calling df.load() on it

        Returns
        -------
        PsiCoverage
            PsiCoverage for prefixes found in all specified files

        Notes
        -----
        Does not check that events are same in each input file. It will fail if
        they are not the same size, which should catch most cases, but be wary
        that events information is derived from the first file alone.
        """

        def backwards_compatible_preprocess(ds: xr.Dataset) -> xr.Dataset:
            """Update PsiCoverage to add missing information from old files"""
            if "prefix_total" not in ds.data_vars.keys():
                get_logger().warning(
                    "%s is an old PsiCoverage file without prefix_total",
                    ds.encoding["source"],
                )
                ds = ds.assign(
                    prefix_total=(
                        "prefix",
                        da.full(
                            ds.sizes["prefix"],
                            np.nan,
                            dtype=ds["raw_total"].dtype,
                            chunks=ds.chunks["prefix"],
                        ),
                    ),
                )
            return ds

        chunks: Dict[str, Optional[int]] = {
            "prefix": -1,
            "ec_idx": -1,
            "bootstrap_replicate": None,
        }
        if ec_idx_nchunks:
            chunks["ec_idx"] = constants.DEFAULT_COVERAGE_CHUNKS * ec_idx_nchunks
        if prefix_nchunks:
            chunks["prefix"] = prefix_nchunks

        if not isinstance(path, list):
            path = [path]
        df = (
            xr.open_mfdataset(
                path,
                engine="zarr",
                group=constants.NC_PSICOVERAGE,
                combine="nested",
                concat_dim="prefix",
                join="override",
                compat="override",
                coords="minimal",
                data_vars="minimal",
                chunks=chunks,
                parallel=True,
                preprocess=backwards_compatible_preprocess,
            )
            .drop_vars("lsv_offsets", errors="ignore")
            .drop_duplicates("prefix")
        )
        if len(path) > 1:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        if preload:
            df.load()
            events_df.load()
        return cls(df, events_df)

    def updated(
        self,
        bootstrap_psi: Optional[xr.DataArray],
        raw_psi: Optional[xr.DataArray],
        **update_attrs,
    ) -> "PsiCoverage":
        """Create updated :py:class:`PsiCoverage` with new values of psi

        Parameters
        ----------
        bootstrap_psi, raw_psi: Optional[xr.DataArray]
            If specified, new values of bootstrap_psi, raw_psi to use in new
            PsiCoverage (with all other variables equal)
        update_attrs:
            Additional kwargs are set as attributes to the dataset used to
            construct the resulting PsiCoverage

        Returns
        -------
        PsiCoverage
            Updated :py:class:`PsiCoverage` with new values of psi
        """
        df = self.df
        # update psi arrays
        if bootstrap_psi is not None:
            if set(self.bootstrap_psi.dims) != set(bootstrap_psi.dims):
                raise ValueError("bootstrap_psi doesn't have same named axes")
            df = df.assign(bootstrap_psi=bootstrap_psi)
        if raw_psi is not None:
            if set(self.raw_psi.dims) != set(raw_psi.dims):
                raise ValueError("raw_psi doesn't have same named axes")
            df = df.assign(raw_psi=raw_psi)
        # update/add attributes
        df = df.assign_attrs(**update_attrs)
        # return resulting PsiCoverage object
        return PsiCoverage(df, self.events_df)

    def _save_df(
        self,
        remove_bam_attrs: bool = False,
    ) -> xr.Dataset:
        """Prepare dataset of psicoverage that will be saved

        Prepare dataset of psicoverage that will be saved. This sets the
        chunking, clears any encodings from before, removes attrs that have to
        do with BAM, etc.
        """
        USE_CHUNKS = dict(
            ec_idx=constants.DEFAULT_COVERAGE_CHUNKS,
            prefix=1,
            bootstrap_replicate=-1,
        )
        save_df = self.df
        if save_df.sizes["ec_idx"] > 0:
            save_df = save_df.chunk(USE_CHUNKS)  # type: ignore[arg-type]
        # remove attributes referring to bam?
        if remove_bam_attrs:
            for k in [x for x in save_df.attrs.keys() if str(x).startswith("bam")]:
                save_df.attrs.pop(k, None)
        return save_df

    def to_zarr(
        self,
        path: Union[str, Path],
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save :py:class:`PsiCoverage` to specified path

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
        save_df = self._save_df()
        save_df_future = cast(
            Delayed,
            save_df.to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICOVERAGE,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = persist_with_progress(save_df_future)
        save_df_future.compute()
        self.events_to_zarr(path, mode="a", consolidated=consolidated)
        return

    def to_zarr_slice(
        self,
        path: Union[str, Path],
        prefix_slice: slice,
    ) -> None:
        """Save :py:class:`PsiCoverage` to specified path for specified slice on prefix

        Save :py:class:`PsiCoverage` to specified path for specified slice.
        Typically run after :py:meth:`PsiCoverage.to_zarr_slice_init`

        Parameters
        ----------
        path: Union[str, Path]
            Path with output file in zarr format with metadata initialized by
            :py:meth:`PsiCoverage.to_zarr_slice_init`
        prefix_slice: slice
            Slice of prefix dimension in output zarr store to save current
            PsiCoverage

        See Also
        --------
        PsiCoverage.to_zarr_slice_init
        """
        self._save_df().drop_vars("prefix").pipe(
            lambda x: x.drop_vars(
                [k for k, v in x.variables.items() if "prefix" not in v.dims]
            )
        ).to_zarr(
            path, group=constants.NC_PSICOVERAGE, region=dict(prefix=prefix_slice)
        )
        return

    @classmethod
    def to_zarr_slice_init(
        cls,
        path: Union[str, Path],
        events_df: xr.Dataset,
        prefixes: Sequence[str],
        num_bootstraps: int,
        cov_dtype: type = np.float32,
        psicov_attrs: Optional[Dict[Hashable, Any]] = None,
    ) -> None:
        """Initialize zarr store for saving :py:class:`PsiCoverage` over many writes

        Initialize zarr for :py:class:`PsiCoverage` over many prefixes.
        Saves all information except dimensions that are prefix-specific.
        This enables multithreaded (or multiprocess) write with
        :py:meth:`PsiCoverage.to_zarr_slice`

        Parameters
        ----------
        path: Union[str, Path]
            Path for output Zarr for psicoverage output
        events_df: xr.Dataset
            Dataset encoding psicoverage events (Events.save_df,
            PsiCoverage.events_df, etc.)
        prefixes: Sequence[str]
            Values for the prefix dimension coordinate
        num_bootstraps: int
            Number of bootstrap replicates that will be used
        cov_dtype: type
            What type to use for psi/total_coverage arrays
        psicov_attrs: Dict[Hashable, Any]
            Attributes to include

        See Also
        --------
        PsiCoverage.to_zarr_slice
        """
        if psicov_attrs is None:
            psicov_attrs = {}
        # force events to be saved as single chunk (no benefit for chunking here)
        events_df = events_df.chunk(events_df.sizes)
        # save events
        events_df.to_zarr(path, mode="w", group=constants.NC_EVENTS, consolidated=False)
        # dims for skeleton
        raw_dims = ("ec_idx", "prefix")
        bootstrap_dims = (*raw_dims, "bootstrap_replicate")
        # shapes for skeleton
        raw_shape = (events_df.sizes["ec_idx"], len(prefixes))
        bootstrap_shape = (*raw_shape, num_bootstraps)
        # chunksizes for skeleton
        raw_chunks = (constants.DEFAULT_COVERAGE_CHUNKS, 1)
        bootstrap_chunks = (*raw_chunks, None)
        # arrays for skeleton
        raw_arr = da.empty(raw_shape, dtype=cov_dtype, chunks=raw_chunks)
        passed_arr = da.empty(raw_shape, dtype=bool, chunks=raw_chunks)
        bootstrap_arr = da.empty(
            bootstrap_shape, dtype=cov_dtype, chunks=bootstrap_chunks
        )
        # save metadata for skeleton
        xr.Dataset(
            dict(
                bootstrap_psi=(bootstrap_dims, bootstrap_arr),
                bootstrap_total=(bootstrap_dims, bootstrap_arr),
                event_passed=(raw_dims, passed_arr),
                raw_psi=(raw_dims, raw_arr),
                raw_total=(raw_dims, raw_arr),
                prefix_total=(
                    "prefix",
                    da.empty(len(prefixes), dtype=cov_dtype, chunks=1),
                ),
            ),
        ).to_zarr(
            path,
            mode="a",
            compute=False,
            group=constants.NC_PSICOVERAGE,
            consolidated=False,
        )
        # save offsets, prefixes, and attributes
        xr.Dataset(
            coords={"prefix": prefixes},
            attrs=psicov_attrs,
        ).to_zarr(path, mode="a", group=constants.NC_PSICOVERAGE, consolidated=True)
        return

    @classmethod
    def convert_sj_batch(
        cls,
        sjs: Sequence[Path],
        lsvs: Events,
        path: Path,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
        imap_unordered_fn=map,
        prefixes: Optional[Sequence[str]] = None,
    ) -> None:
        """Load PsiCoverage from sj paths, save to single output path

        Parameters
        ----------
        sjs: Sequence[Path]
            Paths to input SJ files that will have PsiCoverage evaluated for
        lsvs: Events
            Events over which PsiCoverage will be calculated
        path: Path
            Output path for PsiCoverage zarr file
        minreads, minbins: float
            Quantifiability thresholds
        num_bootstraps: int
            The number of bootstrap replicates for bootstrapped estimates
        pvalue_threshold: float
            P-value threshold for removing stacks under leave-one-out Poisson
            model of per-bin read coverage, for both raw and bootstrapped
            coverage (Set to nonpositive value to skip stack detection)
        imap_unordered_fn
            Loading/saving of input SJ files will be passed through this
            function, which can enable concurrency if appropriate
        prefixes: Optional[Sequence[str]]
            If specified, names for each SJ file that will be stored in the
            resulting file as `prefix`. If None (default), use prefix from BAM
            file name stored in SJ file.
        """
        log = get_logger()

        if prefixes is not None and len(prefixes) != len(sjs):
            raise ValueError(
                "prefixes specified but with different length than sjs"
                f" ({len(sjs) = }, {len(prefixes) = })"
            )

        # how to go from sj path to psicoverage file
        def sj_to_psicov(sj_path: Path) -> PsiCoverage:
            return PsiCoverage.from_sj_lsvs(
                SJExperiment.from_zarr(sj_path),
                lsvs,
                minreads=minreads,
                minbins=minbins,
                num_bootstraps=num_bootstraps,
                pvalue_threshold=pvalue_threshold,
            )

        if len(sjs) == 0:
            raise ValueError("At least one SJ file must be processed")
        elif len(sjs) == 1:
            # if there is only one file, don't bother with imap_unordered_fn
            log.info("Inferring PsiCoverage from %s", sjs[0])
            psi_coverage = sj_to_psicov(sjs[0])
            if prefixes:
                psi_coverage = psi_coverage.rename_prefixes(prefixes)
            log.info("Saving %s to %s", psi_coverage, path)
            psi_coverage.to_zarr(path)
        else:
            # precompute prefixes to use
            log.debug("Precomputing prefixes corresponding to input SJ files")
            prefixes = prefixes or [SJExperiment.prefix_from_zarr(x) for x in sjs]
            # we have more than one input file
            log.info("Saving event information and metadata to %s", path)
            PsiCoverage.to_zarr_slice_init(
                path,
                lsvs.save_df,
                prefixes,
                num_bootstraps,
                psicov_attrs=dict(
                    sj=[str(x) for x in sjs],
                    minreads=minreads,
                    minbins=minbins,
                ),
            )

            def job_fn(sj_idx: int, sj: Path) -> Path:
                sj_to_psicov(sj).to_zarr_slice(
                    path,
                    slice(sj_idx, 1 + sj_idx),
                )
                return sj

            jobs = imap_unordered_fn(lambda x: job_fn(x[0], x[1]), list(enumerate(sjs)))
            for ndx, sj in enumerate(jobs, 1):
                log.info("Saved coverage from %s (%d / %d)", sj, ndx, len(sjs))
        return

    def sum(
        self,
        new_prefix: str,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> "PsiCoverage":
        """Create aggregated :py:class:`PsiCoverage` with sum coverage over prefixes

        Parameters
        ----------
        new_prefix: str
            Prefix for summarized :py:class:`PsiCoverage`
        min_experiments_f: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of prefixes that must pass
            individually for the event to be considered as passed for the group

        Returns
        -------
        PsiCoverage
            Sum coverage over prefixes, with passed being defined over group
            filters
        """
        if self.num_prefixes > 1:
            event_passed = self.passed_min_experiments(min_experiments_f)
            raw_total = self.raw_total.sum("prefix")
            raw_coverage = (self.raw_total * self.raw_psi).sum("prefix")
            raw_psi = (raw_coverage / raw_total.where(raw_total > 0)).fillna(0)
            bootstrap_total = self.bootstrap_total.sum("prefix")
            bootstrap_coverage = (self.bootstrap_total * self.bootstrap_psi).sum(
                "prefix"
            )
            bootstrap_psi = (
                bootstrap_coverage / bootstrap_total.where(bootstrap_total > 0)
            ).fillna(0)
            prefix_total = self.prefix_total.sum("prefix")
            df = xr.Dataset(
                data_vars=dict(
                    event_passed=event_passed,
                    raw_total=raw_total,
                    raw_psi=raw_psi,
                    bootstrap_total=bootstrap_total,
                    bootstrap_psi=bootstrap_psi,
                    prefix_total=prefix_total,
                ),
                attrs=dict(original_prefix=self.prefixes),
            ).expand_dims(prefix=[new_prefix])
        else:
            df = (
                self.df.assign_coords(prefix=[new_prefix])
                    .assign_attrs(original_prefix=self.prefixes)
            )
        return PsiCoverage(df, self.events_df)

    def mask_events(self, passed: xr.DataArray) -> "PsiCoverage":
        """Return :py:class:`PsiCoverage` passing only events that are passed in input

        Parameters
        ----------
        passed: xr.DataArray
            boolean array(ec_idx) where connections marked as not passed
            (False) will be marked as not passed (False) in resulting
            :py:class:`PsiCoverage`

        Returns
        -------
        PsiCoverage
            Same coverage but passing only events that are passed in input (and
            in the original object)
        """
        return self._with_updated_df(
            self.df.assign(event_passed=self.event_passed & passed)
        )

    def dataset(
        self,
        properties: Sequence[str] = constants.DEFAULT_PSI_PROPERTIES,
        quantiles: Sequence[float] = list(),
        psibins: Optional[int] = None,
        any_passed: bool = True,
    ) -> xr.Dataset:
        """Extract selected properties into single :py:class:`xr.Dataset`

        Parameters
        ----------
        properties: Sequence[str]
            PsiCoverage properties to request.
        quantiles: Sequence[float]
            If non-empty, calculate quantiles of posterior distribution
        psibins: Optional[int]
            If specified, calculate discretized approximation to posterior
            distribution with this many bins
        any_passed: bool
            By default, we report any_passed (summarization over event_passed);
            this flag allows this to be excluded.

        Notes
        -----
        quantiles and psibins computations use the approximate posterior
        distribution Beta(approximate_alpha, approximate_beta) because:

        - by default, it's 30 times faster than using the bootstrap mixture,
        - for high coverage events, the bootstrap distribution is discrete (for
          each bootstrap replicate), so the approximate distribution is a
          better representation of our desired model of variability in PSI.
        """
        # initialize variables to return with noting if any experiment passed
        quantify_vars: Dict[str, xr.DataArray] = dict()
        if any_passed:
            quantify_vars["any_passed"] = self.event_passed.any("prefix")
        # add properties
        for x in properties:
            quantify_vars[x] = getattr(self, x)
        if len(quantiles):
            quantify_vars["psi_quantile"] = self.approximate_quantile(quantiles)
        if psibins:
            quantify_vars["psi_pmf"] = self.approximate_discretized_pmf(psibins)
        return xr.Dataset(quantify_vars).reset_coords(drop=True)  # type: ignore[arg-type]

    @classmethod
    def mock_with_psi_and_total(
        cls,
        psi: npt.NDArray[np.float32],
        total: npt.NDArray[np.float32],
        prefix: str = "mock",
        passed: None | npt.NDArray[np.bool_] = None,
        rng: None | np.random.RandomState = None,
    ) -> "PsiCoverage":
        """Returns PsiCoverage over binary events with specified psi/total

        Parameters
        ----------
        psi: ndarray[float]
            1D array with values of raw_psi for first junction in each LSV
        total: ndarray[float]
            1D array with same shape as psi with values of raw_total
        prefix: str
            Prefix associated with input coverage
        passed: Optional[ndarray[bool]]
            1D array with same shape as psi marking if passed. If None, all passed
        rng: Optional[np.random.RandomState]
            random state used to generate bootstrap_psi/total from raw_psi/total.
            If None, uses NumPy global random state (np.random).

        Notes
        -----
        All LSVs are source events with two junctions and no retained intron.
        Bootstrap coverage is determined by Poisson sampling.
        """
        psi = np.asarray(psi)
        total = np.asarray(total)
        if psi.ndim != 1:
            raise ValueError(f"{psi = } must be 1-dimensional")
        if psi.shape != total.shape:
            raise ValueError(f"{psi = } and {total = } must have same shape")
        if (psi < 0).any() or (psi > 1).any():
            raise ValueError(f"{psi = } must have values in [0, 1]")
        if (total < 1).any():
            raise ValueError(f"{total = } must have values >= 1")
        if passed is None:
            passed = np.ones(psi.shape, dtype=np.bool_)
        else:
            passed = np.asarray(passed)
            if psi.shape != passed.shape:
                raise ValueError(
                    f"Not-none value of {passed = } must have same shape as {psi = }"
                )
        # get arrays over ec_idx
        prefix_total = total.astype(np.float32).sum()
        event_passed = np.repeat(passed.astype(np.bool_), 2)
        raw_total = np.repeat(total.astype(np.float32), 2)
        raw_psi = np.empty(raw_total.shape, dtype=np.float32)
        raw_psi[::2] = psi
        raw_psi[1::2] = 1 - psi
        raw_coverage = raw_psi * raw_total
        bootstrap_coverage = (rng.poisson if rng else np.random.poisson)(
            raw_coverage, size=(30, len(raw_psi))
        ).astype(np.float32)
        bootstrap_total = np.repeat(
            np.add.reduceat(bootstrap_coverage, np.arange(0, len(raw_psi), 2), axis=-1),
            2,
            axis=-1,
        )
        bootstrap_psi = bootstrap_coverage / bootstrap_total
        df = (
            xr.Dataset(
                {
                    "bootstrap_psi": (("bootstrap_replicate", "ec_idx"), bootstrap_psi),
                    "bootstrap_total": (
                        ("bootstrap_replicate", "ec_idx"),
                        bootstrap_total,
                    ),
                    "event_passed": ("ec_idx", event_passed),
                    "raw_psi": ("ec_idx", raw_psi),
                    "raw_total": ("ec_idx", raw_total),
                    "prefix_total": prefix_total,
                },
                {},
                {
                    "bam_path": f"mock/{prefix}.bam",
                    "bam_version": nm_version,
                    "minbins": 1.0,
                    "minreads": 1.0,
                },
            )
            .expand_dims({"prefix": [prefix]})
            .transpose("ec_idx", "prefix", "bootstrap_replicate")
        )
        # make events_df
        offsets = np.arange(0, 1 + len(raw_psi), 2).astype(np.uint64)
        connection_idx = np.arange(len(raw_psi)).astype(np.uint64)
        event_type = np.repeat(np.array([b"s"]), len(psi))
        is_intron = np.zeros(len(raw_psi), dtype=np.bool_)
        ref_exon_idx = np.arange(0, len(raw_psi), 2).astype(np.uint64)
        events_df = xr.Dataset(
            {},
            {
                "_offsets": ("e_offsets_idx", offsets),
                "connection_idx": ("ec_idx", connection_idx),
                "event_type": ("e_idx", event_type),
                "is_intron": ("ec_idx", is_intron),
                "ref_exon_idx": ("e_idx", ref_exon_idx),
            },
            {
                "intron_hash": 0,
                "junction_hash": 0,
            },
        )
        return PsiCoverage(df, events_df)
