"""
Exons.py

Exons for genes

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import TYPE_CHECKING, MutableMapping, Optional, Union, cast, List, Sequence

import numpy as np
import numpy.typing as npt
import xarray as xr
import zarr

import rna_majiq.constants as constants
from rna_majiq.internals import AnnotatedTranscripts as _AnnotatedTranscripts

from .GeneRegions import GeneRegions
from .Genes import Genes

from ..logger import get_logger

if TYPE_CHECKING:
    from .GeneIntrons import GeneIntrons
    from .GeneJunctions import GeneJunctions


class AnnotatedTranscripts(GeneRegions):
    """Collection of annotated transcripts per gene and the transcript_id/start/end of each, plus exon start/ends inside each"""

    IDX_NAME = "annotated_transcripts_idx"
    DF_VARS = tuple([*GeneRegions.DF_VARS, "transcript_id", "exon_idx_start", "exon_idx_end"])
    EX_ARRAY_VARS = ("exon_start", "exon_end",)
    ZARR_GROUP = constants.NC_ANNOTATEDTRANSCRIPTS

    def __init__(self, annotatedTranscripts: _AnnotatedTranscripts):
        super().__init__(annotatedTranscripts)
        return

    @property
    def _annotated_transcripts(self) -> _AnnotatedTranscripts:
        """expose underlying internals representation of AnnotatedTranscripts"""
        return self._gene_regions

    @property
    def annotated_transcript_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    @property
    def annotated_transcripts_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    @property
    def transcript_id(self) -> List[str]:
        return self._annotated_transcripts.transcript_id

    @property
    def exon_idx_start(self) -> npt.NDArray[np.int64]:
        return self._annotated_transcripts.exon_idx_start

    @property
    def exon_idx_end(self) -> npt.NDArray[np.int64]:
        return self._annotated_transcripts.exon_idx_end

    @property
    def exon_start(self) -> npt.NDArray[np.int64]:
        return self._annotated_transcripts.exon_start

    @property
    def exon_end(self) -> npt.NDArray[np.int64]:
        return self._annotated_transcripts.exon_end

    @property
    def ex_dfs(self) -> dict:
        out = {}
        for var in self.EX_ARRAY_VARS:
            out[var] = xr.Dataset(coords={var: getattr(self, var)})
        return out

    def to_zarr(
            self,
            store: Union[MutableMapping, str, Path],
            mode: str,
            consolidated: bool = False,
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
        for var, df in self.ex_dfs.items():
            df.to_zarr(
                store,
                mode=mode,
                group=self.ZARR_GROUP,
                consolidated=consolidated
            )

    @classmethod
    def from_arrays(
            cls,
            genes: Genes,
            gene_idx: npt.ArrayLike,
            start: npt.ArrayLike,
            end: npt.ArrayLike,
            transcript_id: npt.ArrayLike,
            exon_idx_start:  npt.ArrayLike,
            exon_idx_endt:  npt.ArrayLike,
            exon_start:  npt.ArrayLike,
            exon_end:  npt.ArrayLike,
    ) -> "AnnotatedTranscripts":
        """Create :class:`AnnotatedTranscripts` from :class:`Genes` and input arrays"""
        start = np.asarray(start)
        end = np.asarray(end)

        return AnnotatedTranscripts(
            _AnnotatedTranscripts(genes._genes, gene_idx, start, end, transcript_id,
                                  exon_idx_start, exon_idx_endt, exon_start, exon_end)
        )

    @classmethod
    def from_zarr(
            cls,
            store: Union[MutableMapping, str, Path],
            genes: Optional[Genes] = None,
    ) -> "AnnotatedTranscripts":
        """Read exons from zarr file/store

        Parameters
        ----------
        store: Union[MutableMapping, str, Path]
            store or path to zarr file
        genes: Optional[Genes]
            genes on which the exons are defined. If None, try loading from
            zarr file. Note that rna_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        if genes is None:
            genes = Genes.from_zarr(store)
        try:
            with xr.open_zarr(store, group=cls.ZARR_GROUP) as df:
                df.load()
                return AnnotatedTranscripts.from_arrays(
                    genes,
                    df.gene_idx.values,
                    df.start.values,
                    df.end.values,
                    df.transcript_id.values.tolist(),
                    df.exon_idx_start.values,
                    df.exon_idx_end.values,
                    df.exon_start.values,
                    df.exon_end.values,
                )
        except zarr.errors.PathNotFoundError:
            log = get_logger()
            log.warning("AnnotatedTranscripts not found in splicegraph, using empty arrays")
            return AnnotatedTranscripts.from_arrays(
                genes,
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
            )


