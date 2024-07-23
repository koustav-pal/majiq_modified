"""
GeneJunctions.py

Container of gene junctions that are open intervals over genes

Author: Joseph K Aicher
"""

from typing import TYPE_CHECKING, Optional

import numpy as np
import numpy.typing as npt

import rna_majiq.constants as constants
from rna_majiq.internals import GeneJunctions as _GeneJunctions

from .Exons import Exons
from .GeneConnections import GeneConnections

if TYPE_CHECKING:
    from .PassedJunctions import GroupJunctionsGenerator, PassedJunctionsGenerator


class GeneJunctions(GeneConnections):
    """Collection of junctions per gene and their coordinates, flags, and exons"""

    IDX_NAME = "gj_idx"
    INTERNALS_CONSTRUCTOR = _GeneJunctions
    ZARR_GROUP = constants.NC_GENEJUNCTIONS

    def __init__(self, gene_junctions: _GeneJunctions):
        super().__init__(gene_junctions)
        return

    def is_denovo(
        self,
        gj_idx: None | int | list[int] | npt.NDArray[np.integer] = None,
        annotated_junctions: Optional["GeneJunctions"] = None,
    ) -> npt.NDArray[np.bool_]:
        """Return denovo status of selected junctions

        Parameters
        ----------
        gj_idx: Optional[array_like[int]]
            Index into junctions to get denovo status for. If None, get denovo
            status for all junctions.
        annotated_junctions: Optional[GeneJunctions]
            If specified, use junctions found in `annotated_junctions` as
            definition of annotated junctions, so that a junction is called
            denovo if and only if it is not found in `annotated_junctions`.
        """
        if not annotated_junctions:
            if gj_idx is None:
                return self.denovo
            else:
                return self.denovo[gj_idx]
        else:
            return ~self.overlaps(annotated_junctions, gj_idx)

    def build_group(self, exons: Exons) -> "GroupJunctionsGenerator":
        """Create :py:class:`GroupJunctionsGenerator` starting from these junctions and exons

        Parameters
        ----------
        exons: Exons
            Exons over the same genes as the junctions, which enable
            identification of most likely genes to which novel junctions belong
        """
        from .PassedJunctions import GroupJunctionsGenerator

        return GroupJunctionsGenerator(self, exons)

    def builder(self) -> "PassedJunctionsGenerator":
        """Create :py:class:`PassedJunctionsGenerator` starting from these junctions"""
        from .PassedJunctions import PassedJunctionsGenerator

        return PassedJunctionsGenerator(self)

    @property
    def _gene_junctions(self) -> _GeneJunctions:
        """Underlying internals representation"""
        return self._gene_connections

    @property
    def gj_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx
