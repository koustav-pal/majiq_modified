"""
GeneConnections.py

Parent class that wraps rna_majiq internals for gene connections (introns,
junctions)

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import (
    Callable,
    ClassVar,
    MutableMapping,
    Optional,
    Sequence,
    Type,
    TypeVar,
    Union,
)

import numpy as np
import numpy.typing as npt
import xarray as xr

from .Exons import Exons, _Exons
from .GeneRegions import GeneRegions
from .Genes import Genes

SelfT = TypeVar("SelfT", bound="GeneConnections")


class GeneConnections(GeneRegions):
    DF_VARS = tuple(
        [
            *GeneRegions.DF_VARS,
            "denovo",
            "passed_build",
            "simplified",
            "start_exon_idx",
            "end_exon_idx",
        ]
    )
    ZARR_GROUP = "connections"
    # used by from_arrays()
    INTERNALS_CONSTRUCTOR: ClassVar[Callable] = lambda *args, **kwargs: None  # noqa: E731

    def __init__(self, gene_connections):
        super().__init__(gene_connections)
        return

    @property
    def _gene_connections(self):
        """Underlying internals class for gene connections"""
        return self._gene_regions

    def checksum(self):
        """Checksum including passed/simplified/connections status"""
        return self._gene_connections.checksum()

    def checksum_nodata(self):
        """Checksum only considering gene_idx, start, and end"""
        return self._gene_connections.checksum_nodata()

    def connect_exons(self, exons: Exons) -> None:
        """Connect regions to specified exons"""
        self._gene_connections.connect_exons(exons._exons)
        return

    @property
    def connected_exons(self) -> Optional[Exons]:
        """exons the connections are associated with (or None otherwise)"""
        raw: Optional[_Exons] = self._gene_connections.connected_exons
        return None if raw is None else Exons(raw)

    @property
    def denovo(self) -> npt.NDArray[np.bool_]:
        """Indicate if each connection is denovo or not"""
        return self._gene_connections.denovo

    @property
    def passed_build(self) -> npt.NDArray[np.bool_]:
        """Indicate if each connection passed build filters (reliable) or not"""
        return self._gene_connections.passed_build

    @property
    def simplified(self) -> npt.NDArray[np.bool_]:
        """Indicate if each connection is simplified or not"""
        return self._gene_connections.simplified

    @property
    def start_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Indicate exon_idx associated with start coordinate"""
        return self._gene_connections.start_exon_idx

    @property
    def end_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Indicate exon_idx associated with end coordinate"""
        return self._gene_connections.end_exon_idx

    def src_exon_idx(
        self, region_idx: Optional[npt.ArrayLike] = None
    ) -> npt.NDArray[np.uint64]:
        if region_idx is None:
            region_idx = self._region_idx
        return self._gene_connections.src_exon_idx(region_idx)

    def dst_exon_idx(
        self, region_idx: Optional[npt.ArrayLike] = None
    ) -> npt.NDArray[np.uint64]:
        if region_idx is None:
            region_idx = self._region_idx
        return self._gene_connections.dst_exon_idx(region_idx)

    def _pass_all(self) -> None:
        """Set all connections to have passed build"""
        self._gene_connections._pass_all()
        return

    def _simplify_all(self) -> None:
        """Set all connections to the simplified state"""
        self._gene_connections._simplify_all()
        return

    def _unsimplify_all(self) -> None:
        """Set all connections to the unsimplified state"""
        self._gene_connections._unsimplify_all()
        return

    def to_zarr(
        self,
        store: Union[MutableMapping, str, Path],
        mode: str,
        consolidated: bool = True,
        drop_vars: Optional[Union[str, Sequence[str]]] = None,
    ) -> None:
        """Serialize to zarr format. Note parent genes need to be saved separately.

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            Save to this path/store using Zarr group `self.ZARR_GROUP`
        mode: str
            Indicate mode for writing data to file
        consolidated: bool
            If True, have xarray/zarr consolidate metadata for the entire Zarr
            store.
        drop_vars: Optional[Union[str, Sequence[str]]]
            Drop these variables from `self.df` when writing to file.
            If None (default), drop the index variable, start_exon_idx, and
            end_exon_idx.
        """
        drop_vars = drop_vars or [self.IDX_NAME, "start_exon_idx", "end_exon_idx"]
        super().to_zarr(store, mode, consolidated=consolidated, drop_vars=drop_vars)
        return

    @classmethod
    def from_zarr(
        cls: Type[SelfT],
        store: Union[MutableMapping, str, Path],
        genes: Optional[Genes] = None,
    ) -> SelfT:
        """Load connections from zarr file or store.

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            store or path to zarr file
        genes: Optional[Genes]
            genes on which the connections are defined.
            If None, try loading from zarr store.
            Note that rna_majiq checks if objects refer to the same genes (not
            that they are identical), so it is usually desired to provide the
            variable than using the default behavior.
        """
        if genes is None:
            genes = Genes.from_zarr(store)
        with xr.open_zarr(store, group=cls.ZARR_GROUP) as df:
            df.load()
            return cls.from_arrays(
                genes,
                df.gene_idx.values,
                df.start.values,
                df.end.values,
                denovo=df.denovo.values,
                passed_build=df.passed_build.values,
                simplified=df.simplified.values,
            )

    @classmethod
    def from_arrays(
        cls: Type[SelfT],
        genes: Genes,
        gene_idx: npt.ArrayLike,
        start: npt.ArrayLike,
        end: npt.ArrayLike,
        denovo: Optional[npt.ArrayLike] = None,
        passed_build: Optional[npt.ArrayLike] = None,
        simplified: Optional[npt.ArrayLike] = None,
        connected_exons: Optional[Exons] = None,
    ) -> SelfT:
        """Create connections object from :class:`Genes` and input arrays"""
        if denovo is None:
            denovo = np.zeros_like(gene_idx, dtype=np.bool_)
        if passed_build is None:
            passed_build = np.zeros_like(gene_idx, dtype=np.bool_)
        if simplified is None:
            simplified = np.zeros_like(gene_idx, dtype=np.bool_)
        return cls(
            cls.INTERNALS_CONSTRUCTOR(
                genes._genes,
                gene_idx,
                start,
                end,
                denovo,
                passed_build,
                simplified,
                connected_exons=getattr(connected_exons, "_exons", None),
            )
        )
