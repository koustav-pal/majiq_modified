"""
SJIntrons.py

contig regions not intersected by gene exons (stranded or unstranded), noting
if intersecting annotated intron or not (part of annotated exon)

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import MutableMapping, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import rna_majiq.constants as constants
from rna_majiq.internals import SJIntrons as _SJIntrons

from .ContigRegions import ContigRegions
from .Contigs import Contigs
from .Exons import Exons
from .GeneIntrons import GeneIntrons


class SJIntrons(ContigRegions):
    IDX_NAME = "si_idx"
    DF_VARS = tuple([*ContigRegions.DF_VARS, "annotated"])
    ZARR_GROUP = constants.NC_SJINTRONS

    def __init__(self, sj_introns: _SJIntrons):
        super().__init__(sj_introns)
        return

    @property
    def _sj_introns(self) -> _SJIntrons:
        """expose underlying internals representation of SJIntrons"""
        return self._contig_regions

    @property
    def si_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    @property
    def annotated(self) -> npt.NDArray[np.bool_]:
        """Indicates if region associated with annotated intron"""
        return self._sj_introns.annotated

    @classmethod
    def from_exons_and_introns(
        cls, exons: Exons, gene_introns: GeneIntrons, stranded: bool
    ) -> "SJIntrons":
        """Construct SJIntrons from exons and gene introns"""
        return SJIntrons(
            _SJIntrons.from_exons_and_introns(
                exons._exons, gene_introns._gene_introns, stranded
            )
        )

    @classmethod
    def from_arrays(
        cls,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
        annotated: npt.ArrayLike,
    ) -> "SJIntrons":
        """Create :class:`SJIntrons` from :class:`Contigs` and input arrays"""
        return SJIntrons(
            _SJIntrons(contigs._contigs, contig_idx, start, end, strand, annotated)
        )

    @classmethod
    def from_contigs(cls, contigs: Contigs) -> "SJIntrons":
        """Empty SJIntrons linked to specified contigs"""
        return SJIntrons.from_arrays(
            contigs,
            np.array([], dtype=np.uint64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype="S1"),
            np.array([], dtype=bool),
        )

    @classmethod
    def from_zarr(
        cls,
        store: Union[MutableMapping, str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "SJIntrons":
        """Read SJIntrons from zarr file/store

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            store or path to zarr file
        contigs: Optional[Contigs]
            contigs on which the introns are defined. If None, try loading from
            zarr file.
        """
        if contigs is None:
            contigs = Contigs.from_zarr(store)
        with xr.open_zarr(store, group=cls.ZARR_GROUP) as df:
            df.load()
            return SJIntrons.from_arrays(
                contigs,
                df.contig_idx.values,
                df.start.values,
                df.end.values,
                df.strand.values,
                df.annotated.values,
            )
