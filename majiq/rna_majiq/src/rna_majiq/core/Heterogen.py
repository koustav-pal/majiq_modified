"""
Heterogen.py

Quantify/test differences between two groups of independent experiments

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import TYPE_CHECKING, Collection, Final, Optional, Sequence, Union

import numpy as np
import xarray as xr

import rna_majiq.constants as constants
from rna_majiq.logger import get_logger

from ..core.MixinHasEvents import MixinHasEvents
from .PsiCoverage import PsiCoverage
from .PsiGroup import PsiGroup

if TYPE_CHECKING:
    from ..voila.HeterogenDataset import HeterogenDataset


class Heterogen(MixinHasEvents):
    """Compare Psi between two groups of PsiCoverage (independence assumption)

    Compare Psi between two groups of PsiCoverage under the assumption that
    each experiment is independent. Perform null-hypothesis testing on
    posterior means and/or samples under assumption that underlying values of
    PSI per experiment are independent and identically distributed (i.i.d.).

    Parameters
    ----------
    psi1, psi2: Union[PsiCoverage, PsiGroup]
        :class:`PsiCoverage` or :class:`PsiGroup` for the two groups to compare
    min_experiments_f: float
        Number or proportion of experiments required to pass in each group for
        quantification
    name1, name2: str
        Names to indicate group identity
    show_progress: bool
        Show progress of intermediate computations required for setting up this
        object.
    """

    def __init__(
        self,
        psi1: Union[PsiCoverage, PsiGroup],
        psi2: Union[PsiCoverage, PsiGroup],
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        name1: str = "grp1",
        name2: str = "grp2",
        show_progress: bool = False,
    ):
        log = get_logger()
        # check that name1 != name2
        if name1 == name2:
            raise ValueError(f"{name1 = } and {name2 = } must be different")
        # check events matching
        if psi1.num_connections != psi2.num_connections:
            raise ValueError(
                "psi1 and psi2 must have the same number of connections"
                f" ({psi1.num_connections=}, {psi2.num_connections=})"
            )
        if psi1.num_events != psi2.num_events:
            raise ValueError(
                "psi1 and psi2 must have the same number of events"
                f" ({psi1.num_events=}, {psi2.num_events=})"
            )
        # save events
        MixinHasEvents.__init__(self, None, psi1.events_df)
        # check prefixes overlap
        if prefix_overlap := set(psi1.prefixes) & set(psi2.prefixes):
            # warning if there is overlap between prefixes
            log.warning(
                f"Heterogen input groups have overlapping prefixes ({prefix_overlap})"
            )
        # get events that passed both psi1 and psi2 in enough experiments
        passed = psi1.passed_min_experiments(
            min_experiments_f
        ) & psi2.passed_min_experiments(min_experiments_f)
        self.psi1_original_prefix = psi1.prefixes
        self.psi2_original_prefix = psi2.prefixes
        psi1 = psi1.mask_events(passed)
        psi2 = psi2.mask_events(passed)
        if not isinstance(psi1, PsiGroup):
            log.info("Creating temporary intermediates for group %s", name1)
            psi1 = psi1.group(show_progress=show_progress)
        if not isinstance(psi2, PsiGroup):
            log.info("Creating temporary intermediates for group %s", name2)
            psi2 = psi2.group(show_progress=show_progress)
        self.psi1: Final[PsiGroup] = psi1
        self.psi2: Final[PsiGroup] = psi2
        self.passed: Final[xr.DataArray] = passed
        self.name1: Final[str] = name1
        self.name2: Final[str] = name2
        return

    @cached_property
    def psi_concatenated(self) -> PsiGroup:
        """Concatenated values from `self.psi1` and `self.psi2` for statistics"""
        # since name1 and name2 are unique, will not overlap
        result = PsiGroup.concat(
            self.psi1.rename_prefixes(
                [f"{self.name1}-{x}" for x in self.psi1.prefixes]
            ),
            self.psi2.rename_prefixes(
                [f"{self.name2}-{x}" for x in self.psi2.prefixes]
            ),
        )
        return result

    @property
    def labels(self) -> xr.DataArray:
        """labels associated with prefixes in `self.psi_concatenated`"""
        return xr.DataArray(
            np.repeat([True, False], [self.psi1.num_prefixes, self.psi2.num_prefixes]),
            [("prefix", self.psi_concatenated.prefixes)],
        )

    def raw_stats(
        self,
        ec_idx: Optional[slice] = None,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    ) -> xr.Dataset:
        """Statistics on means, samples from raw posteriors"""
        indexer_kwargs = dict(ec_idx=ec_idx) if ec_idx else dict()
        return self.psi_concatenated.raw_stats(
            self.labels, use_stats=use_stats, **indexer_kwargs
        )

    def approximate_stats(
        self,
        ec_idx: Optional[slice] = None,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        drop: bool = True,
    ) -> xr.Dataset:
        """Statistics on means, samples from approximate posteriors"""
        indexer_kwargs = dict(ec_idx=ec_idx) if ec_idx else dict()
        return self.psi_concatenated.approximate_stats(
            self.labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=drop,
            **indexer_kwargs,
        )

    def dataset(
        self,
        pvalue_quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
    ) -> "HeterogenDataset":
        from ..voila.HeterogenDataset import HeterogenDataset

        return HeterogenDataset.from_heterogen(
            self,
            pvalue_quantiles=pvalue_quantiles,
            use_stats=use_stats,
            psisamples=psisamples,
        )
