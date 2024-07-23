"""
DeltaPsi.py

Quantify differences between two groups of replicate experiments

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import TYPE_CHECKING, Final, Optional, cast

import numpy as np
import xarray as xr

import rna_majiq.beta_mixture as bm
import rna_majiq.constants as constants

from .DPsiPrior import DPsiPrior
from .MixinPsiInference import _logsumexp
from .PMFSummaries import PMFSummaries
from .PsiCoverage import PsiCoverage

if TYPE_CHECKING:
    from ..voila.DeltaPsiDataset import DeltaPsiDataset


class DeltaPsiPMF(PMFSummaries):
    """Specialization of PMFSummaries for DeltaPsi on [-1, 1]

    Parameters
    ----------
    p: xr.DataArray
        Probability distribution(s) over intervals (dimension: pmf_bin)
        Dimensions: [..., "pmf_bin"]
        Coordinates: pmf_bin_start, pmf_bin_end (..., "pmf_bin")
        First pmf_bin_start must be -1, Last pmf_bin_end must be 1
    """

    def __init__(self, p: xr.DataArray):
        """Initialize :class:`DeltaPsiPMF` with given probability vector

        Parameters
        ----------
        p: xr.DataArray
            Probability distribution(s) over intervals (dimension: pmf_bin)
            Dimensions: [..., "pmf_bin"]
            Coordinates: pmf_bin_start, pmf_bin_end (..., "pmf_bin")
            First pmf_bin_start must be -1, Last pmf_bin_end must be 1
        """
        PMFSummaries.__init__(self, p)
        if self.bin_start.values[0] != -1:
            raise ValueError("First bin must start at -1")
        if self.bin_end.values[-1] != 1:
            raise ValueError("Last bin must end at 1")
        return

    def probability_nonchanging(
        self,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
    ) -> xr.DataArray:
        """Probability that abs(dPSI) <= nonchanging_threshold"""
        return self.interval_probability(-nonchanging_threshold, nonchanging_threshold)

    def probability_changing(
        self,
        changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD,
    ) -> xr.DataArray:
        """Probability that abs(dPSI) > changing_threshold"""
        return 1 - self.interval_probability(-changing_threshold, changing_threshold)


class DeltaPsi(object):
    """Compute DeltaPsi between two groups of PsiCoverage (replicate assumption)

    Compute DeltaPsi between two groups of PsiCoverage under assumption that
    the experiments in each group are replicates with the given DPsiPrior.

    Parameters
    ----------
    psi1, psi2: PsiCoverage
        PsiCoverage files to be treated as replicates
    prior: DPsiPrior
        Prior on deltapsi for inference
    min_experiments_f: float
        Only pass events that pass this many experiments in both groups
    name1, name2: str
        Names to use for the aggregated groups coverage/PSI
    """

    def __init__(
        self,
        psi1: PsiCoverage,
        psi2: PsiCoverage,
        prior: DPsiPrior,
        psibins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        name1: str = "psi1",
        name2: str = "psi2",
    ):
        if psi1.num_connections != psi2.num_connections:
            raise ValueError(
                "psi1 and psi2 must have the same number of connections"
                f" ({psi1.num_connections=}, {psi2.num_connections=})"
            )
        # save aggregate psi coverage, original prior
        self.psi1: Final[PsiCoverage] = psi1.sum(
            name1, min_experiments_f=min_experiments_f
        )
        self.psi2: Final[PsiCoverage] = psi2.sum(
            name2, min_experiments_f=min_experiments_f
        )
        self.prior: Final[DPsiPrior] = prior
        self.psibins: Final[int] = psibins
        return

    @property
    def num_connections(self) -> int:
        return self.psi1.num_connections

    @property
    def name1(self) -> str:
        return self.psi1.prefixes[0]

    @property
    def name2(self) -> str:
        return self.psi2.prefixes[0]

    def __repr__(self) -> str:
        return (
            f"DeltaPsi[{self.num_connections}] for {self.name1} vs {self.name2}"
            f" computed with {self.psibins} psi bins and {self.prior}"
        )

    def rebin(self, psibins: int) -> "DeltaPsi":
        """Get DeltaPsi with different bins"""
        return DeltaPsi(
            self.psi1,
            self.psi2,
            self.prior,
            psibins=psibins,
            name1=self.name1,
            name2=self.name2,
        )

    @cached_property
    def passed(self) -> xr.DataArray:
        return self.psi1.event_passed.squeeze(
            "prefix", drop=True
        ) & self.psi2.event_passed.squeeze("prefix", drop=True)

    @cached_property
    def discrete_logprior(self) -> xr.DataArray:
        return self.prior.discretized_logpmf(psibins=self.psibins).astype(
            self.psi1.bootstrap_psi.dtype
        )

    def smooth_logposterior(
        self, ec_idx: Optional[slice] = None, nchunks: int = 128
    ) -> xr.DataArray:
        """Infer deltapsi logposterior on approximation of psi bootstrapped mixture

        Parameters
        ----------
        ec_idx: Optional[slice]
            If specified, only compute logposterior on given slice
        nchunks: int
            Number of chunks to split computation over. The computation is
            expensive, so this allows Dask to split the computation over
            independent threads/workers and give progress updates as each chunk
            completes.
            The result will not have chunks over ec_idx.
        """
        ds = xr.Dataset(
            dict(
                a1=self.psi1.approximate_alpha.squeeze("prefix", drop=True).expand_dims(
                    mix1=1
                ),
                b1=self.psi1.approximate_beta.squeeze("prefix", drop=True).expand_dims(
                    mix1=1
                ),
                a2=self.psi2.approximate_alpha.squeeze("prefix", drop=True).expand_dims(
                    mix2=1
                ),
                b2=self.psi2.approximate_beta.squeeze("prefix", drop=True).expand_dims(
                    mix2=1
                ),
            )
        ).where(self.passed)
        if ec_idx:
            ds = ds.isel(ec_idx=ec_idx)
        nchunks = max(1, nchunks)
        # divide number of event connections into this many chunks for computation
        ds = ds.chunk({"ec_idx": -(ds.sizes["ec_idx"] // -nchunks)})
        return xr.apply_ufunc(
            bm.dpsi_discrete,
            ds["a1"],
            ds["b1"],
            ds["a2"],
            ds["b2"],
            self.discrete_logprior,
            input_core_dims=[["mix1"], ["mix1"], ["mix2"], ["mix2"], ["pmf_bin"]],
            output_core_dims=[["pmf_bin"]],
            dask="allowed",
        ).chunk({"ec_idx": -1})

    def bootstrap_logposterior(
        self, ec_idx: Optional[slice] = None, nchunks: int = 128
    ) -> xr.DataArray:
        """Logposterior average of deltapsi inference over each bootstrap replicate

        Parameters
        ----------
        ec_idx: Optional[slice]
            If specified, only compute logposterior on given slice
        nchunks: int
            Number of chunks to split computation over. The computation is
            expensive, so this allows Dask to split the computation over
            independent threads/workers and give progress updates as each chunk
            completes.
            The result will not have chunks over ec_idx.
        """
        ds = xr.Dataset(
            dict(
                a1=self.psi1.bootstrap_alpha.rename(prefix="mix1"),
                b1=self.psi1.bootstrap_beta.rename(prefix="mix1"),
                a2=self.psi2.bootstrap_alpha.rename(prefix="mix2"),
                b2=self.psi2.bootstrap_beta.rename(prefix="mix2"),
            )
        ).where(self.passed)
        if ec_idx:
            ds = ds.isel(ec_idx=ec_idx)
        nchunks = max(1, nchunks)
        # divide number of event connections into this many chunks for computation
        ds = ds.chunk({"ec_idx": -(ds.sizes["ec_idx"] // -nchunks)})
        per_bootstrap = xr.apply_ufunc(
            bm.dpsi_discrete,
            ds["a1"],
            ds["b1"],
            ds["a2"],
            ds["b2"],
            self.discrete_logprior,
            input_core_dims=[["mix1"], ["mix1"], ["mix2"], ["mix2"], ["pmf_bin"]],
            output_core_dims=[["pmf_bin"]],
            dask="allowed",
        )
        # take sum over bootstrap replicates, which isn't normalized
        if per_bootstrap.size:
            unnormalized = xr.apply_ufunc(
                _logsumexp,
                per_bootstrap,
                input_core_dims=[["bootstrap_replicate"]],
                dask="allowed",
            )
        else:
            # reduction over per_bootstrap with 0-dimension can lead to error
            # in reductionwith logsumexp, so use a simpler reduction that
            # doesn't fail
            unnormalized = per_bootstrap.sum(dim=["bootstrap_replicate"])
        # normalize each resulting distribution
        if unnormalized.size:
            logZ = xr.apply_ufunc(
                _logsumexp,
                unnormalized,
                input_core_dims=[["pmf_bin"]],
                dask="allowed",
            )
        else:
            # simpler reduction for zero-length array
            logZ = unnormalized.sum(dim=["pmf_bin"])
        return (unnormalized - logZ).chunk({"ec_idx": -1})

    @cached_property
    def discrete_prior(self) -> DeltaPsiPMF:
        """:class:`DeltaPsiPMF` for discretized deltapsi prior"""
        return DeltaPsiPMF(cast(xr.DataArray, np.exp(self.discrete_logprior)))

    @cached_property
    def smooth_posterior(self) -> DeltaPsiPMF:
        """:class:`DeltaPsiPMF` for dpsi posterior from approximated PSI posteriors"""
        return DeltaPsiPMF(cast(xr.DataArray, np.exp(self.smooth_logposterior())))

    @cached_property
    def bootstrap_posterior(self) -> DeltaPsiPMF:
        """:class:`DeltaPsiPMF` for average bootstrapped dpsi posteriors"""
        return DeltaPsiPMF(cast(xr.DataArray, np.exp(self.bootstrap_logposterior())))

    def dataset(self, nchunks: int = 128) -> "DeltaPsiDataset":
        """Reduce to :class:`DeltaPsiDataset` used for VOILA visualization"""
        from ..voila.DeltaPsiDataset import DeltaPsiDataset

        return DeltaPsiDataset.from_deltapsi(self, nchunks=nchunks)
