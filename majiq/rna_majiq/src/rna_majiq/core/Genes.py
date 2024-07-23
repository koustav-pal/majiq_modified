"""
Genes.py

Container of genes that are closed intervals over different contigs

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import List, MutableMapping, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import rna_majiq.constants as constants
from rna_majiq.internals import Genes as _Genes

from .ContigRegions import ContigRegions
from .Contigs import Contigs


class Genes(ContigRegions):
    """Collection of genes on :class:`Contigs` with their coordinates, and ids

    Parameters
    ----------
    genes: _Genes
        Underlying object binding the internal C++ API

    See Also
    --------
    Contigs.from_zarr
    """

    IDX_NAME = "gene_idx"
    DF_VARS = tuple([*ContigRegions.DF_VARS, "gene_id", "gene_name"])
    ZARR_GROUP = constants.NC_GENES

    def __init__(self, genes: _Genes):
        super().__init__(genes)
        return

    def checksum(self):
        return self._genes.checksum()

    @property
    def _genes(self) -> _Genes:
        """expose underlying internals representation of Genes"""
        return self._contig_regions

    @property
    def gene_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    def __getitem__(self, gene_id: str) -> int:
        """Get gene_idx for specified gene_id"""
        try:
            return self._genes[gene_id]
        except IndexError as err:
            raise KeyError(f"{gene_id = } not found in Genes instance") from err

    @property
    def gene_id(self) -> List[str]:
        return self._genes.gene_id

    @property
    def gene_name(self) -> List[str]:
        return self._genes.gene_name

    def annotate_gene_idx(self, df: xr.Dataset) -> xr.Dataset:
        """For now, just add gene_id to df using df.gene_idx"""
        return df.assign_coords(
            gene_id=(df.gene_idx.dims, np.array(self.gene_id)[df.gene_idx]),
        )

    @classmethod
    def from_arrays(
        cls,
        contigs: Contigs,
        contig_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        strand: npt.ArrayLike,
        gene_id: List[str],
        gene_name: Optional[List[str]] = None,
    ) -> "Genes":
        """Create :class:`Genes` from :class:`Contigs` and input arrays

        Parameters
        ----------
        contigs: Contigs
            Contigs on which genes are defined
        contig_idx: array_like, int
            Index of contig for each gene
        start, end: array_like, int
            Coordinates for each gene on the contig
        strand: array_like, {"+", "-"}
            Forward or reverse strand for each gene
        gene_id: List[str]
            Unique identifier for each gene
        gene_name: Optional[List[str]]
            Name for each gene. If not specified, will be the same as gene_id

        Notes
        -----
        Genes must be provided in sorted order (contig_idx, start, end, strand).
        """
        if gene_name is None:
            gene_name = gene_id
        return Genes(
            _Genes(contigs._contigs, contig_idx, start, end, strand, gene_id, gene_name)
        )

    @classmethod
    def from_zarr(
        cls,
        store: Union[MutableMapping, str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "Genes":
        """Load :class:`Genes` from specified path/store

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            Store or path to zarr file
        contigs: Optional[Contigs]
            contigs on which the genes are defined. If None, try loading from
            zarr file. Note that rna_majiq checks if objects refer to the
            same contigs (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        if contigs is None:
            contigs = Contigs.from_zarr(store)
        with xr.open_zarr(store, group=cls.ZARR_GROUP) as df:
            df.load()
            return Genes.from_arrays(
                contigs,
                df.contig_idx.values,
                df.start.values,
                df.end.values,
                df.strand.values,
                df.gene_id.values.tolist(),
                df.gene_name.values.tolist(),
            )
