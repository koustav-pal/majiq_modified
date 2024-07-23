"""
SJJunctions.py

contig regions not intersected by gene exons (stranded or unstranded), noting
if intersecting annotated junction or not (part of annotated exon)

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import MutableMapping, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import rna_majiq.constants as constants
from rna_majiq.internals import SJJunctions as _SJJunctions

from .ContigRegions import ContigRegions
from .Contigs import Contigs


class SJJunctions(ContigRegions):
    IDX_NAME = "sj_idx"
    ZARR_GROUP = constants.NC_SJJUNCTIONS

    def __init__(self, sj_junctions: _SJJunctions):
        super().__init__(sj_junctions)
        return

    @property
    def _sj_junctions(self) -> _SJJunctions:
        """expose underlying internals representation of SJJunctions"""
        return self._contig_regions

    def to_unstranded(self) -> "SJJunctions":
        """Get unique junctions ignoring strand (mark all unstranded)"""
        return SJJunctions(self._sj_junctions.to_unstranded())

    def flip_strand(self) -> "SJJunctions":
        """Get junctions with strands flipped in sorted order"""
        return SJJunctions(self._sj_junctions.flip_strand())

    @property
    def sj_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    @classmethod
    def from_arrays(
        cls,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
    ) -> "SJJunctions":
        """Create :class:`SJJunctions` from :class:`Contigs` and input arrays"""
        return SJJunctions(
            _SJJunctions(contigs._contigs, contig_idx, start, end, strand)
        )

    @classmethod
    def from_zarr(
        cls,
        store: Union[MutableMapping, str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "SJJunctions":
        """Read SJJunctions from zarr file/store

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            store or path to zarr file
        contigs: Optional[Contigs]
            contigs on which the junctions are defined. If None, try loading
            from zarr file.
        """
        if contigs is None:
            contigs = Contigs.from_zarr(store, group=constants.NC_SJJUNCTIONSCONTIGS)
        with xr.open_zarr(store, group=cls.ZARR_GROUP) as df:
            df.load()
            return SJJunctions.from_arrays(
                contigs,
                df.contig_idx.values,
                df.start.values,
                df.end.values,
                df.strand.values,
            )
