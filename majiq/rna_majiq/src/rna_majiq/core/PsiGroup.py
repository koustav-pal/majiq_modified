"""
PsiGroup.py

Rechunked attributes from :class:`PsiCoverage` so that prefixes are in a single
chunk, connections are chunked

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import (
    Dict,
    Final,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
    overload,
)

import xarray as xr
import zarr
from dask.delayed import Delayed
from moccasin.moccasin import persist_with_progress

import rna_majiq.constants as constants

from .PsiCoverage import MixinPsiOverPrefixes, PsiCoverage


class PsiGroup(MixinPsiOverPrefixes):
    """PSI posterior distributions parameters for a group of experiments.

    Raw and approximate bootstrap posterior distributions on PSI for a group of
    experiments, chunked appropriately for summaries, operations, etc. over the
    groups.
    """

    EXPECTED_VARIABLES = {
        "raw_alpha": ("ec_idx", "prefix"),
        "raw_beta": ("ec_idx", "prefix"),
        "approximate_alpha": ("ec_idx", "prefix"),
        "approximate_beta": ("ec_idx", "prefix"),
        "prefix": ("prefix",),
    }

    def __init__(
        self,
        df: xr.Dataset,
        events: xr.Dataset,
        hold_temporary: Tuple[TemporaryDirectory, ...] = tuple(),
    ):
        # save PSI information and events
        MixinPsiOverPrefixes.__init__(self, df, events)
        # hold onto temporary directory (if any)
        self.hold_temporary: Final = hold_temporary
        return

    @classmethod
    def concat(
        cls,
        *objs,
        override_args: Optional[Sequence] = None,
        update_kwargs: Optional[Dict] = None,
    ):
        """Concatenate multiple :class:`PsiGroup` into single `PsiGroup`

        Parameters
        ----------
        *objs:
            Combine the data from these input objects, using additional
            args/kwargs from the first element. At least one item must be
            specified.
        override_args: Optional[Sequence]
            Instead of using `self._non_df_args` for additional positional
            arguments, use these positional arguments
        update_kwargs: Dict
            Add additional keyword arguments to constructor.
            If hold_temporary is specified, will concatenate temporary
            directories in addition to those found in `objs`.

        Notes
        -----
        Overrides functionality from :class:`MixinSubsettablePrefixes` in order
        to ensure that all parent temporary directories are held
        """
        if update_kwargs is None:
            update_kwargs = {}
        update_kwargs = update_kwargs.copy()
        update_kwargs["hold_temporary"] = tuple(
            [x for obj in objs for x in obj.hold_temporary]
            + list(update_kwargs.get("hold_temporary", tuple()))
        )
        return super().concat(
            *objs, override_args=override_args, update_kwargs=update_kwargs
        )

    @property
    def _non_df_kwargs(self):
        return dict(hold_temporary=self.hold_temporary)

    def mask_events(self, passed: xr.DataArray) -> "PsiGroup":
        """Return :class:`PsiGroup` passing only events passed in input

        Parameters
        ----------
        passed: xr.DataArray
            boolean array where unpassed connections will be masked out over
            all prefixes

        Returns
        -------
        PsiGroup
            Posterior distributions for connections where passed=True (NaN
            otherwise)
        """
        return self._with_updated_df(self.df.where(passed))

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Save :class:`PsiGroup` to specified path as zarr"""
        self.df.chunk(
            {"ec_idx": constants.DEFAULT_PSIGROUP_CHUNKS, "prefix": -1}
        ).to_zarr(path, mode="w", group=constants.NC_PSIGROUP, consolidated=False)
        self.events_to_zarr(path, mode="a", consolidated=True)
        return

    # overload from_zarr_permissive to indicate that if always_convert is given
    # a literal of True, result will always be PsiGroup
    @overload
    @classmethod
    def from_zarr_permissive(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        always_convert: Literal[True] = ...,
        chunk_like_group: bool = ...,
        ec_idx_nchunks: Optional[int] = ...,
        show_progress: bool = ...,
    ) -> "PsiGroup":
        ...

    @overload
    @classmethod
    def from_zarr_permissive(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        always_convert: Literal[False] = ...,
        chunk_like_group: bool = ...,
        ec_idx_nchunks: Optional[int] = ...,
        show_progress: bool = ...,
    ) -> Union["PsiGroup", PsiCoverage]:
        ...

    @overload
    @classmethod
    def from_zarr_permissive(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        always_convert: bool = ...,
        chunk_like_group: bool = ...,
        ec_idx_nchunks: Optional[int] = ...,
        show_progress: bool = ...,
    ) -> Union["PsiGroup", PsiCoverage]:
        ...

    @classmethod
    def from_zarr_permissive(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        always_convert: bool = False,
        chunk_like_group: bool = True,
        ec_idx_nchunks: Optional[int] = 1,
        show_progress: bool = False,
    ):
        """Load :class:`PsiGroup` or :class:`PsiCoverage` from Zarr paths.

        Load :class:`PsiGroup` or :class:`PsiCoverage` from Zarr paths.
        Will convert to `PsiGroup` unless all input paths are for
        :class:`PsiCoverage` and `always_convert` is False.

        Parameters
        ----------
        path: Union[str, Path, List[Union[str, Path]]]
            Path(s) with group information saved in zarr format.
            Can be PsiCoverage or PsiGroup saved in the zarr store.
        always_convert: bool
            If True, the result will always be `PsiGroup`.
        chunk_like_group: bool
            If the result is `PsiCoverage`, should it (True) be opened to chunk
            like `PsiGroup` (to use as is for computation), or (False) should
            it be opened so as to efficiently be converted later?
        ec_idx_nchunks: Optional[int]
            Number of chunks to open along connections relative to PsiGroup
            chunksize.
        show_progress: bool
            Show progress if needing to convert :class:`PsiCoverage` to
            `PsiGroup`.
        """
        # make path a list
        if not isinstance(path, list):
            path = [path]
        # identify valid PsiGroup files:
        events_groups_match = [
            (constants.NC_EVENTS in x, constants.NC_PSIGROUP in x)
            for x in (zarr.open(p, mode="r") for p in path)
        ]
        # they must all have events in them
        if missing := [
            p
            for p, (has_events, has_group) in zip(path, events_groups_match)
            if not has_events
        ]:
            raise ValueError(
                f"Input paths {missing} appear to be missing event information"
            )
        has_group = [has_group for has_events, has_group in events_groups_match]

        # get paths with groups information
        group_path = [p for p, x in zip(path, has_group) if x]
        # assume that if it doesn't have PsiGroup, that it has PsiCoverage
        psicov_path = [p for p, x in zip(path, has_group) if not x]

        if not psicov_path:
            return PsiGroup.from_zarr(group_path, ec_idx_nchunks=ec_idx_nchunks)
        elif not (group_path or always_convert):
            prefix_nchunks: Optional[int]
            if chunk_like_group:
                if ec_idx_nchunks is not None and ec_idx_nchunks > 0:
                    ec_idx_nchunks = max(
                        1,
                        (ec_idx_nchunks * constants.DEFAULT_PSIGROUP_CHUNKS)
                        // constants.DEFAULT_COVERAGE_CHUNKS,
                    )
                prefix_nchunks = None
            else:
                ec_idx_nchunks = None
                prefix_nchunks = 1
            return PsiCoverage.from_zarr(
                psicov_path,
                ec_idx_nchunks=ec_idx_nchunks,
                prefix_nchunks=prefix_nchunks,
            )
        # otherwise, we definitely have psicov that needs to be converted
        psicov_group = PsiCoverage.from_zarr(
            psicov_path, ec_idx_nchunks=None, prefix_nchunks=1
        ).group(ec_idx_nchunks=ec_idx_nchunks, show_progress=show_progress)
        if group_path:
            return PsiGroup.concat(
                psicov_group,
                PsiGroup.from_zarr(group_path, ec_idx_nchunks=ec_idx_nchunks),
            )
        else:
            return psicov_group

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        ec_idx_nchunks: Optional[int] = 1,
        hold_temporary: Optional[TemporaryDirectory] = None,
    ) -> "PsiGroup":
        """Load :class:`PsiGroup` from specified Zarr path.

        Parameters
        ----------
        path: Union[str, Path, List[Union[str, Path]]]
            Path(s) with group information saved in zarr format
        ec_idx_nchunks: Optional[int]
            Number of chunks over ec_idx to load simultaneously at once.
            If None, load all chunks along the dimension together.
        hold_temporary: Optional[TemporaryDirectory]
            Keep a temporary directory alive (typically used from
            :meth:`PsiGroup.from_psicov`)
        """
        chunks: Dict[str, Optional[int]] = {
            "ec_idx": -1,
            "prefix": -1,
        }
        if ec_idx_nchunks:
            chunks["ec_idx"] = ec_idx_nchunks * constants.DEFAULT_PSIGROUP_CHUNKS
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_PSIGROUP,
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="no_conflicts",
            coords="minimal",
            data_vars="minimal",
            chunks=chunks,
            parallel=True,
        ).drop_duplicates("prefix")
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        return PsiGroup(
            df,
            events_df,
            hold_temporary=tuple() if hold_temporary is None else (hold_temporary,),
        )

    @classmethod
    def from_psicov(
        cls,
        psicov: PsiCoverage,
        save_zarr: Optional[Union[str, Path]] = None,
        tmp_base_dir: Optional[Union[str, Path]] = None,
        ec_idx_nchunks: Optional[int] = 1,
        show_progress: bool = False,
    ) -> "PsiGroup":
        """Create :class:`PsiGroup` from :class:`PsiCoverage`.

        Create class with `EXPECTED_VARIABLES` from :class:`PsiCoverage` with
        prefixes in single chunk and connections chunked as appropriate saved
        to permanent or temporary storage.

        Parameters
        ----------
        psicov: PsiCoverage
            Coverage for group of experiments.
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
        """
        ec_chunksize = constants.DEFAULT_PSIGROUP_CHUNKS
        hold_temporary: Optional[TemporaryDirectory] = None
        if save_zarr is None:
            hold_temporary = TemporaryDirectory(dir=tmp_base_dir)
            save_zarr = Path(hold_temporary.name) / "zarr"

        with TemporaryDirectory(dir=tmp_base_dir) as intermediate_dir:
            intermediate_path = Path(intermediate_dir) / "zarr"
            intermediate_ds = xr.Dataset(
                {
                    "raw_alpha": psicov.raw_alpha,
                    "raw_beta": psicov.raw_beta,
                    "approximate_alpha": psicov.approximate_alpha,
                    "approximate_beta": psicov.approximate_beta,
                },
                {
                    "prefix": psicov.prefixes,
                },
            )
            intermediate_future = cast(
                Delayed,
                intermediate_ds.to_zarr(intermediate_path, compute=False),
            )
            # save to intermediate file (always temporary)
            if show_progress:
                intermediate_future = persist_with_progress(intermediate_future)
            intermediate_future.compute()
            # load back into memory, rechunk, and save to save_zarr
            params_ds = xr.open_zarr(intermediate_path).chunk(
                {"ec_idx": ec_chunksize, "prefix": -1}
            )
            for v in params_ds.variables.values():
                v.encoding.clear()
            params_future = cast(
                Delayed,
                params_ds.to_zarr(
                    save_zarr,
                    mode="w",
                    group=constants.NC_PSIGROUP,
                    compute=False,
                    consolidated=False,
                ),
            )
            if show_progress:
                params_future = persist_with_progress(params_future)
            params_future.compute()
        psicov.events_to_zarr(save_zarr, mode="a", consolidated=True)
        return PsiGroup.from_zarr(
            save_zarr, ec_idx_nchunks=ec_idx_nchunks, hold_temporary=hold_temporary
        )

    @property
    def raw_alpha(self) -> xr.DataArray:
        return self.df["raw_alpha"]

    @property
    def raw_beta(self) -> xr.DataArray:
        return self.df["raw_beta"]

    @property
    def approximate_alpha(self) -> xr.DataArray:
        return self.df["approximate_alpha"]

    @property
    def approximate_beta(self) -> xr.DataArray:
        return self.df["approximate_beta"]

    @cached_property
    def event_passed(self) -> xr.DataArray:
        """array(ec_idx, prefix) indicating if event passed"""
        return self.raw_alpha.notnull()
