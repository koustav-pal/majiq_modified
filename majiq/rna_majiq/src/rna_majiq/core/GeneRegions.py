"""
ContigRegions.py

Parent class that wraps rna_majiq internals for gene regions

Author: Joseph K Aicher
"""

import numpy as np
import numpy.typing as npt

from .Genes import Genes
from .Regions import Regions


class GeneRegions(Regions):
    DF_VARS = tuple(["gene_idx", *Regions.DF_VARS])

    def __init__(self, gene_regions):
        super().__init__(gene_regions)
        return

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}["
            f"{len(self.genes.contigs)} contigs, {len(self.genes)} genes, {len(self)}]"
        )

    @property
    def _gene_regions(self):
        """Underlying internals class for gene regions"""
        return self._regions

    @property
    def genes(self) -> Genes:
        """Genes for which these regions are defined"""
        return Genes(self._parents)

    @property
    def gene_idx(self) -> npt.NDArray[np.uint64]:
        """Index of gene on which region is defined"""
        return self._gene_regions.gene_idx

    def index(
        self,
        gene_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
    ) -> npt.NDArray[np.int64]:
        """Get indexes of gene regions (-1 if not present)"""
        return self._gene_regions.index(gene_idx, start, end)

    def slice_for_gene(self, gene_idx: int) -> slice:
        return self._slice_for_parent(gene_idx)
