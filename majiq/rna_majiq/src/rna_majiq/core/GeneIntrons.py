"""
GeneIntrons.py

Container of gene introns that are open intervals over genes

Author: Joseph K Aicher
"""

from typing import TYPE_CHECKING, Optional

import numpy as np
import numpy.typing as npt

import rna_majiq.constants as constants
from rna_majiq.internals import GeneIntrons as _GeneIntrons

from .Exons import Exons
from .GeneConnections import GeneConnections
from .Genes import Genes

if TYPE_CHECKING:
    from .GroupIntronsGenerator import GroupIntronsGenerator


class GeneIntrons(GeneConnections):
    """Collection of introns per gene and their coordinates, flags, and exons

    Parameters
    ----------
    gene_introns: _GeneIntrons
        Underlying object binding the internal C++ API
    """

    IDX_NAME = "gi_idx"
    INTERNALS_CONSTRUCTOR = _GeneIntrons
    ZARR_GROUP = constants.NC_GENEINTRONS

    def __init__(self, gene_introns: _GeneIntrons):
        super().__init__(gene_introns)
        return

    def propagate_to_annotated(
        self,
        annotated_exons: Optional[Exons] = None,
        keep_annotated: bool = constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        discard_denovo: bool = not constants.DEFAULT_BUILD_DENOVO_IR,
    ) -> "GeneIntrons":
        """Get :class:`GeneIntrons` for annotated exons

        Parameters
        ----------
        annotated_exons: Optional[Exons]
            Propagate introns to match annotated exons. If not specified, get
            annotated exons from connected exons (raise error if not connected
            exons).
        keep_annotated: bool
            Keep all annotated introns regardless of whether they passed
        discard_denovo: bool
            Discard all denovo introns regardless of whether they passed
        """
        if not annotated_exons:
            if self.connected_exons:
                annotated_exons = self.connected_exons.get_annotated()
            else:
                raise ValueError(
                    "annotated_exons must be passed when introns not connected"
                )
        return (
            annotated_exons.potential_introns(make_simplified=True)
            .update_flags_from(self)
            .filter_passed(keep_annotated=keep_annotated, discard_denovo=discard_denovo)
        )

    def propagate_through_annotated(
        self,
        keep_annotated: bool = constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        discard_denovo: bool = not constants.DEFAULT_BUILD_DENOVO_IR,
    ) -> "GeneIntrons":
        """Propagate introns to annotated exons, then back to current exons.

        This ensures that all possible introns between annotated exons have the
        same flags, which keeps introns consistent with two-pass build.
        This requires that the introns are connected to exons.
        """
        if not self.connected_exons:
            raise ValueError("GeneIntrons are not connected to exons")
        return self.propagate_to_annotated(
            keep_annotated=True, discard_denovo=False
        ).propagate_to_annotated(
            annotated_exons=self.connected_exons,
            keep_annotated=keep_annotated,
            discard_denovo=discard_denovo,
        )

    def is_denovo(
        self,
        gi_idx: None | int | list[int] | npt.NDArray[np.integer] = None,
        annotated_introns: Optional["GeneIntrons"] = None,
    ) -> npt.NDArray[np.bool_]:
        """Return denovo status of selected introns

        Parameters
        ----------
        gi_idx: Optional[array_like[int]]
            Index into introns for to get denovo status for.
            If None, get denovo status for all introns.
        annotated_introns: Optional[GeneIntrons]
            If specified, use introns that are found in `annotated_exons` to
            reduce the number of exons called annotated by treating exons
            overlapping with `annotated_exons` as annotated.
            If `annotated_introns` has connected exons, will appropriately
            propagate to original gaps between annotated exons.
        """
        if not annotated_introns:
            if gi_idx is None:
                return self.denovo
            else:
                return self.denovo[gi_idx]
        else:
            # propagate annotated_introns to potential introns between
            # underlying annotated exons (if possible)
            try:
                annotated_introns = annotated_introns.propagate_to_annotated()
            except ValueError:
                pass  # keep as original value
            return ~self.overlaps(annotated_introns, gi_idx)

    def build_group(self) -> "GroupIntronsGenerator":
        """Create :py:class:`GroupIntronsGenerator` to update these introns in place"""
        from .GroupIntronsGenerator import GroupIntronsGenerator

        return GroupIntronsGenerator(self)

    def filter_passed(
        self,
        keep_annotated: bool = constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        discard_denovo: bool = not constants.DEFAULT_BUILD_DENOVO_IR,
    ) -> "GeneIntrons":
        """Return :py:class:`GeneIntrons` subset that all passed build filters

        Parameters
        ----------
        keep_annotated: bool
            Keep all annotated introns regardless of whether they passed
        discard_denovo: bool
            Discard all denovo introns regardless of whether they passed
        """
        return GeneIntrons(
            self._gene_introns.filter_passed(keep_annotated, discard_denovo)
        )

    def update_flags_from(self, donor_introns: "GeneIntrons") -> "GeneIntrons":
        """Update flags using overlapping donor :py:class:`GeneIntrons`

        Update flags (denovo, passed_build, simplified) using overlapping donor
        :py:class:`GeneIntrons` in place.

        Parameters
        ----------
        donor_introns: GeneIntrons
            Flags from donor_introns will be used to update flags in self when
            they overlap (if annotated in donor, mark as annotated, if passed
            in donor, mark as passed, if unsimplified in donor, mark as
            unsimplified)

        Returns
        -------
        GeneIntrons
            Returns self, after flags have been updated
        """
        self._gene_introns.update_flags_from(donor_introns._gene_introns)
        return self

    @property
    def _gene_introns(self) -> _GeneIntrons:
        """Underlying internals representation"""
        return self._gene_connections

    @property
    def gi_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    @classmethod
    def from_genes(
        cls, genes: Genes, connected_exons: Optional[Exons] = None
    ) -> "GeneIntrons":
        """Empty introns matched to specified genes"""
        return GeneIntrons.from_arrays(
            genes,
            np.array([], dtype=np.uint64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            denovo=np.array([], dtype=bool),
            passed_build=np.array([], dtype=bool),
            simplified=np.array([], dtype=bool),
            connected_exons=connected_exons,
        )
