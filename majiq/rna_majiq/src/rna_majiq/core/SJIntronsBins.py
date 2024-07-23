"""
SJIntronsBins.py

Bins with raw read coverage from input experiments

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import MutableMapping, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import rna_majiq.constants as constants
from rna_majiq._version import version as nm_version
from rna_majiq.internals import ExperimentStrandness
from rna_majiq.internals import SJIntronsBins as _SJIntronsBins
from rna_majiq.logger import get_logger

from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .SJBinsReads import SJBinsReads
from .SJIntrons import SJIntrons


class SJIntronsBins(SJBinsReads):
    """Per-bin read coverage over introns

    Parameters
    ----------
    sj_intronsbins: _SJIntronsBins
        Underlying object binding the internal C++ API
    strandness: ExperimentStrandness
        strandness that was used for the experiment
    original_path: str
        Original path to BAM file from which coverage was derived
    original_version, original_time: str
        Version of MAJIQ and time BAM file was processed
    gene_introns_checksum, exons_checksum: Optional[int]
        Checksum of gene introns, exons used to define regions for coverage

    Notes
    -----
    Intron coverage is dependent on the gene introns and exons used to define
    regions on which to assess coverage.

    See Also
    --------
    SJExperiment.from_bam
    SJIntronsBins.from_bam
    SJIntronsBins.from_zarr
    """

    IDX_NAME = "sib_idx"

    def __init__(
        self,
        sj_intronsbins: _SJIntronsBins,
        strandness: ExperimentStrandness,
        original_path: str,
        original_version: str,
        original_time: str,
        prefix: Optional[str] = None,
        gene_introns_checksum: Optional[int] = None,
        exons_checksum: Optional[int] = None,
        **additional_metadata,
    ):
        super().__init__(
            sj_intronsbins,
            strandness,
            original_path,
            original_version,
            original_time,
            prefix=prefix,
            gene_introns_checksum=gene_introns_checksum,
            exons_checksum=exons_checksum,
            **additional_metadata,
        )
        return

    @property
    def gene_introns_checksum(self) -> Optional[int]:
        """checksum of introns used in constructing :class:`SJIntronsBins`"""
        return self._additional_metadata["gene_introns_checksum"]

    @property
    def exons_checksum(self) -> Optional[int]:
        """checksum of introns used in constructing :class:`SJIntronsBins`"""
        return self._additional_metadata["exons_checksum"]

    @property
    def _sj_intronsbins(self) -> _SJIntronsBins:
        return self._sj_binsreads

    @property
    def sib_idx(self):
        return self._sjbin_idx

    @property
    def regions(self) -> SJIntrons:
        return SJIntrons(self._regions)

    @property
    def sib_idx_start(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_start

    @property
    def sib_idx_end(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_end

    @classmethod
    def from_bam(
        cls,
        path: Union[str, Path],
        total_bins: int,
        exons: Exons,
        gene_introns: GeneIntrons,
        strandness: ExperimentStrandness = constants.DEFAULT_BAM_STRANDNESS,
        nthreads: int = constants.DEFAULT_BAM_NTHREADS,
    ) -> "SJIntronsBins":
        """Load SJIntronsBins from BAM file

        Parameters
        ----------
        path: Union[str, Path]
            Path to input BAM file
        total_bins: int
            Maximum possible bin for intron coverage (usually obtained from
            SJJunctionsBins)
        exons: Exons
            exons used to define sj_introns on which to infer coverage
        gene_introns: GeneIntrons
            introns used to identify sj_introns that can only count for
            annotated introns
        strandness: ExperimentStrandness
            The strand-specificity of the experiment
        nthreads: int
            Number of threads to use to read in BAM file
        """
        path = str(Path(path).resolve())
        original_version = nm_version
        original_time = str(np.datetime64("now"))
        return SJIntronsBins(
            _SJIntronsBins.from_bam(
                path,
                total_bins,
                exons._exons,
                gene_introns._gene_introns,
                strandness,
                nthreads,
            ),
            strandness,
            path,
            original_version,
            original_time,
            gene_introns_checksum=gene_introns.checksum_nodata(),
            exons_checksum=exons.checksum(),
        )

    def to_zarr(
        self,
        store: Union[MutableMapping, str, Path],
        consolidated: bool = True,
        check_experiment_if_exists: bool = True,
    ) -> None:
        """Serialize to zarr format

        Will only write to existing file if has SJJunctionsBins with same
        original path and version

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            Store or Path for output zarr file
        consolidated: bool
            Should the zarr file be consolidated (e.g. we don't expect array
            metadata to change) at the end of function?
        check_experiment_if_exists: bool
            If zarr file already exists at path, check that it has
            SJJunctionsBins and that the original path and version match
        """
        if check_experiment_if_exists:
            try:
                check_experiment = Path(store).exists()  # type: ignore[arg-type]
            except TypeError:
                # store is not pathlike, so assume to be store
                check_experiment = True
        else:
            check_experiment = False
        if check_experiment:
            try:
                with xr.open_zarr(store, group=constants.NC_SJJUNCTIONSBINS) as x:
                    if x.original_path != self.original_path:
                        raise ValueError(
                            f"Will not save SJIntronsBins to existing file {store}"
                            " which appears to have junctions reads from"
                            f" path {x.original_path}"
                            f" vs introns from path {self.original_path}"
                        )
                    if x.original_version != self.original_version:
                        raise ValueError(
                            f"Will not save SJIntronsBins to existing file {store}"
                            " which appears to have junctions reads from"
                            f" version {x.original_version}"
                            f" vs introns from version {self.original_version}"
                        )
            except (KeyError, OSError, AttributeError) as err:
                raise ValueError(
                    f"Will not save SJIntronsBins to existing file {store}"
                    " since it does not appear to have matching junctions reads"
                ) from err
            # otherwise, appears to have compatible SJJunctionsBins which
            # implies that they have matching contigs
        # otherwise
        self.regions.contigs.to_zarr(store, "w", consolidated=False)
        self.regions.to_zarr(store, "a", consolidated=False)
        self._df.pipe(lambda x: x.chunk(x.sizes)).to_zarr(
            store,
            mode="a",
            group=constants.NC_SJINTRONSBINS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_arrays(
        cls,
        regions: SJIntrons,
        bin_reads: npt.ArrayLike,
        bin_idx: npt.ArrayLike,
        offsets: npt.ArrayLike,
        total_bins: int,
        strandness: ExperimentStrandness,
        original_path: str = "<none>",
        original_version: str = nm_version,
        original_time: Optional[str] = None,
        prefix: Optional[str] = None,
        gene_introns_checksum: Optional[int] = None,
        exons_checksum: Optional[int] = None,
    ) -> "SJIntronsBins":
        """Create :class:`SJIntronsBins` from :class:`SJIntrons`, arrays, metadata"""
        if original_time is None:
            original_time = str(np.datetime64("now"))
        return SJIntronsBins(
            _SJIntronsBins(
                regions._sj_introns, bin_reads, bin_idx, offsets, total_bins
            ),
            strandness,
            original_path,
            original_version,
            original_time,
            prefix=prefix,
            gene_introns_checksum=gene_introns_checksum,
            exons_checksum=exons_checksum,
        )

    @classmethod
    def from_regions(
        cls,
        introns: SJIntrons,
        total_bins: int,
        strandness: ExperimentStrandness = ExperimentStrandness.NONE,
        original_path: str = "<none>",
        original_version: str = nm_version,
        original_time: Optional[str] = None,
    ) -> "SJIntronsBins":
        """Empty SJIntronsBins matched to input introns"""
        return SJIntronsBins.from_arrays(
            introns,
            np.array([], dtype=np.float32),
            np.array([], dtype=np.int32),
            np.zeros(1 + len(introns), dtype=np.uint64),
            total_bins,
            strandness,
            original_path=original_path,
            original_version=original_version,
            original_time=original_time,
        )

    @classmethod
    def from_zarr(cls, store: Union[MutableMapping, str, Path]) -> "SJIntronsBins":
        """Load SJIntronsBins from zarr format"""
        regions = SJIntrons.from_zarr(store)
        with xr.open_zarr(store, group=constants.NC_SJINTRONSBINS) as df:
            df.load()
            try:
                strandness = ExperimentStrandness(ord(df.strandness[0]))
            except AttributeError:
                get_logger().warning(
                    f"SJJunctionsBins in {store} did not save strandness"
                    " -> defaulting to NONE"
                )
                strandness = ExperimentStrandness.NONE
            return SJIntronsBins.from_arrays(
                regions,
                df.bin_reads.values,
                df.bin_idx.values,
                df._offsets.values,
                df.total_bins,
                strandness,
                original_path=df.original_path,
                original_version=df.original_version,
                original_time=df.original_time,
                prefix=df.attrs.get("prefix"),
                gene_introns_checksum=df.attrs.get("gene_introns_checksum"),
                exons_checksum=df.attrs.get("exons_checksum"),
            )
