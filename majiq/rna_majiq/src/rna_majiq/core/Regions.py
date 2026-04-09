"""
Regions.py

Parent class that wraps rna_majiq internals for regions. This is a parent class

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import ClassVar, Final, MutableMapping, Optional, Sequence, Union

import numpy as np
import numpy.typing as npt
import xarray as xr


class Regions(object):
    IDX_NAME: ClassVar[str] = "_region_idx"
    DF_VARS: ClassVar[Sequence[str]] = ("start", "end")
    ZARR_GROUP: ClassVar[str] = "regions"

    def _get_df(self, idx_name: str, use_vars: Sequence[str]) -> xr.Dataset:
        return xr.Dataset(coords={v: (idx_name, getattr(self, v)) for v in use_vars})

    @property
    def df(self) -> xr.Dataset:
        """xr.Dataset of specified data"""
        return self._get_df(self.IDX_NAME, [self.IDX_NAME, *self.DF_VARS])

    def __init__(self, regions):
        self._regions: Final = regions
        return

    def to_zarr(
        self,
        store: Union[MutableMapping, str, Path],
        mode: str,
        consolidated: bool = True,
        drop_vars: Optional[Union[str, Sequence[str]]] = None,
    ) -> None:
        """Serialize to zarr format. Note parent regions need to be saved separately

        Parameters
        ----------
        path: Union[MutableMapping, str, Path]
            Save to this store/path using Zarr group `self.ZARR_GROUP`
        mode: str
            Indicate mode for writing data to store
        consolidated: bool
            If True, have xarray/zarr consolidate metadata for the entire Zarr
            store.
            This increases speed of subsequent loads, but should only be done
            when saving last group to the zarr store.
            That is, if writing contigs, then genes, then exons to the same
            path, set conslidated=True only on the write for exons.
        drop_vars: Optional[Union[str, Sequence[str]]]
            Drop these variables from `self.df` when writing to file.
            If None (default), drop only the index variable.
            Will raise error if the variables do not exist.
            Will not complain if it makes it no longer possible to reload the
            file because a necessary variable was dropped -- primarily for
            internal use.
        """
        self.df.drop_vars(drop_vars or self.IDX_NAME).pipe(
            lambda x: x.chunk(x.sizes)
        ).to_zarr(
            store,
            mode=mode,
            group=self.ZARR_GROUP,
            consolidated=consolidated,
        )

    def __eq__(self, other):
        try:
            return self._regions == other._regions
        except AttributeError:
            return False

    def __len__(self) -> int:
        """Number of regions"""
        return len(self._regions)

    def overlaps(
        self, other, region_idx: Optional[npt.ArrayLike] = None
    ) -> npt.NDArray[np.bool_]:
        """Get mask over region_idx indicating if they overlap regions in other

        Parameters
        ----------
        other:
            Regions of same type as self
        region_idx: Optional[array_like[int]]
            Indexes for regions in self. If None, compute overlaps for all
            regions in self

        Returns
        -------
        array[bool]
            array matching region_idx. An element is true when there exists a
            feature in other that overlaps self[region_idx].

        Notes
        -----
        For junctions, this just searches for a matching feature.
        Requires that the regions share the same parent.
        """
        if region_idx is None:
            region_idx = np.arange(len(self))
        return self._regions.overlaps(region_idx, other._regions)

    @property
    def _region_idx(self) -> npt.NDArray[np.int64]:
        """Index over regions"""
        return np.arange(len(self))

    @property
    def _parents(self):
        """internals class for parents on which regions defined"""
        return self._regions._parents

    @property
    def _parent_idx_start(self) -> npt.NDArray[np.uint64]:
        """First index in regions for each parent"""
        return self._regions._parent_idx_start

    @property
    def _parent_idx_end(self) -> npt.NDArray[np.uint64]:
        """One after last index in regions for each parent"""
        return self._regions._parent_idx_end

    def _slice_for_parent(self, parent_idx: int) -> slice:
        """Get slice into regions for specified parent"""
        import sys; print(f"DEBUG parent_idx={parent_idx} start={self._parent_idx_start[parent_idx]} end={self._parent_idx_end[parent_idx]}", file=sys.stderr, flush=True)
        return slice(
            self._parent_idx_start[parent_idx],
            self._parent_idx_end[parent_idx],
        )

    @property
    def start(self) -> npt.NDArray[np.int64]:
        """Start coordinate of each region"""
        return self._regions.start

    @property
    def end(self) -> npt.NDArray[np.int64]:
        """End coordinate of each region"""
        return self._regions.end
