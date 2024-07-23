"""
MixinPsiInference.py

Abstract mixin classes to define how different PSI quantities are computed
from alpha, beta, etc.

Author: Joseph K Aicher
"""

from abc import ABC, abstractmethod
from functools import cached_property
from typing import Collection, Mapping, Optional, Sequence, Tuple, TypeVar, Union, cast

import numpy as np
import xarray as xr

import rna_majiq.beta_mixture as bm
import rna_majiq.constants as constants
from rna_majiq._stats import nanmedian, nanquantile

from ..experiments import min_experiments
from .MixinHasEvents import MixinHasEvents
from .MixinSubsettablePrefixes import MixinSubsettablePrefixes


class MixinRawPsi(ABC):
    """Methods for PSI inference when raw_alpha, raw_beta are defined"""

    @property
    @abstractmethod
    def raw_alpha(self) -> xr.DataArray:
        """array(...) of alpha for raw posterior"""
        ...

    @property
    @abstractmethod
    def raw_beta(self) -> xr.DataArray:
        """array(...) of beta for raw posterior"""
        ...

    @cached_property
    def raw_alpha_plus_beta(self) -> xr.DataArray:
        return self.raw_alpha + self.raw_beta

    @cached_property
    def raw_total(self) -> xr.DataArray:
        """Raw total for event assuming generalized Jeffreys prior.

        Raw total for event assuming generalized Jeffreys prior.
        MAJIQ uses a prior that always adds to 1.
        Therefore, given posterior parameters, we can infer the original raw
        total by taking their sum and subtracting 1 (the sum of the prior
        parameters).
        """
        return self.raw_alpha_plus_beta - 1

    @cached_property
    def raw_posterior_mean(self) -> xr.DataArray:
        """array(...) means of raw posterior distribution on PSI"""
        return self.raw_alpha / (self.raw_alpha_plus_beta)

    @cached_property
    def raw_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of raw posterior distribution on PSI"""
        mean = self.raw_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.raw_alpha_plus_beta)

    @cached_property
    def raw_posterior_std(self) -> xr.DataArray:
        """array(...) standard deviations of raw posterior distribution"""
        return cast(xr.DataArray, np.sqrt(self.raw_posterior_variance))

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        """array(...) means of raw posterior distribution on PSI (alias)"""
        return self.raw_posterior_mean

    @property
    def raw_psi_variance(self) -> xr.DataArray:
        """array(...) variances of raw posterior distribution on PSI (alias)"""
        return self.raw_posterior_variance

    @property
    def raw_psi_std(self) -> xr.DataArray:
        """array(...) standard deviations of raw posterior distribution (alias)"""
        return self.raw_posterior_std

    @cached_property
    def _raw_alpha_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.raw_alpha
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def _raw_beta_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.raw_beta
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    def raw_stats(
        self,
        labels: xr.DataArray,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        **indexer_kwargs: slice,
    ) -> xr.Dataset:
        """Statistics on raw posterior means with respect to `labels`

        Parameters
        ----------
        labels: xr.DataArray
            Boolean labels over prefixes to perform tests on.
            Note that if it shares any dimensions with optional
            `indexer_kwargs`, those slices will be applied to labels as well as
            posterior distribution parameters
        use_stats: Union[str, Collection[str]]
            Names of test statistics to perform.
            Must be in `constants.STATS_AVAILABLE`.
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case testing is done
            over the remaining prefixes.
            All input dimensions are expected in distribution parameters.
            Slice will be applied to `labels` whenever the dimension is found
            in `labels`.

        Returns
        -------
        xr.Dataset
            With variable `pvalue` for statistics on raw posterior means
        """
        quantiles: Sequence[float] = []
        psisamples: int = 0
        drop: bool = True
        a = self._raw_alpha_core_prefix
        b = self._raw_beta_core_prefix
        if indexer_kwargs:
            a = a.isel(indexer_kwargs)
            b = b.isel(indexer_kwargs)
            labels = labels.isel(indexer_kwargs, missing_dims="ignore")
        return _compute_stats(
            a,
            b,
            labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=drop,
        )


class MixinApproximatePsi(ABC):
    """Methods for PSI inference when approximate_alpha, approximate_beta defined"""

    @property
    @abstractmethod
    def approximate_alpha(self) -> xr.DataArray:
        """array(...) of alpha for smooth posterior (approximation of bootstrap mixture)"""
        ...

    @property
    @abstractmethod
    def approximate_beta(self) -> xr.DataArray:
        """array(...) of beta for smooth posterior (approximation of bootstrap mixture)"""
        ...

    @cached_property
    def approximate_alpha_plus_beta(self) -> xr.DataArray:
        return self.approximate_alpha + self.approximate_beta

    @cached_property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        """array(...) means of bootstrap posterior distribution on PSI"""
        return self.approximate_alpha / self.approximate_alpha_plus_beta

    @cached_property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of bootstrap posterior distribution on PSI"""
        mean = self.bootstrap_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.approximate_alpha_plus_beta)

    @cached_property
    def bootstrap_posterior_std(self) -> xr.DataArray:
        """array(...) standard deviations of bootstrap posterior distribution"""
        return cast(xr.DataArray, np.sqrt(self.bootstrap_posterior_variance))

    @property
    def bootstrap_psi_mean(self) -> xr.DataArray:
        """array(...) means of bootstrap posterior distribution on PSI (alias)"""
        return self.bootstrap_posterior_mean

    @property
    def bootstrap_psi_variance(self) -> xr.DataArray:
        """array(...) variances of bootstrap posterior distribution on PSI (alias)"""
        return self.bootstrap_posterior_variance

    @property
    def bootstrap_psi_std(self) -> xr.DataArray:
        """array(...) standard deviations of bootstrap posterior distribution (alias)"""
        return self.bootstrap_posterior_std

    def approximate_cdf(
        self,
        x: Union[xr.DataArray, Sequence[float]],
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute cdf of approximate/smoothed bootstrapped posterior

        Parameters
        ----------
        x: Union[xr.DataArray, Sequence[float]]
            potential realizations of distribution to evaluate CDF at
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects

        Returns
        -------
        xr.DataArray
            Return array(...) of cdf probabilities per connection. If
            `x` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "x"
        """
        alpha = self.approximate_alpha
        beta = self.approximate_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_cdf(alpha, beta, x)

    def approximate_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = (0.1, 0.9),
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute quantiles of approximate/smoothed bootstrapped posterior

        Parameters
        ----------
        quantiles: Union[xr.DataArray, Sequence[float]]
            quantiles of distribution to compute
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects

        Returns
        -------
        xr.DataArray
            Return array(...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"
        """
        alpha = self.approximate_alpha
        beta = self.approximate_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_quantile(alpha, beta, quantiles=quantiles)

    def approximate_discretized_pmf(
        self,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        midpoint_approximation: bool = True,
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute discretized PMF of approximate/smoothed bootstrap posterior

        Parameters
        ----------
        nbins: int
            Number of uniform bins on [0, 1] on which probability mass will be
            computed
        midpoint_approximation: bool
            If midpoint_approximation, approximate PMF by computing
            unnormalized PDF over midpoints of uniform bins, then normalizing.
            Otherwise, use incomplete beta function.
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects
        """
        alpha = self.approximate_alpha
        beta = self.approximate_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_discretized_pmf(
            alpha, beta, midpoint_approximation=midpoint_approximation, nbins=nbins
        )

    @cached_property
    def _approximate_alpha_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.approximate_alpha
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def _approximate_beta_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.approximate_beta
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    def approximate_stats(
        self,
        labels: xr.DataArray,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        drop: bool = True,
        **indexer_kwargs: slice,
    ) -> xr.Dataset:
        """Statistics on approximate posterior means and psisamples.

        Statistics on approximation of bootstrap posterior means and quantiles
        of statistics from repeated sampling (psisamples) with respect to
        `labels`.

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles of pvalues to extract from repeating test statistics on
            `psisamples` draws from the posterior distributions
        psisamples: int
            Number of repeated draws from the posterior distributions to
            perform test statistics on
        labels: xr.DataArray
            Boolean labels over prefixes to perform tests on.
            Note that if it shares any dimensions with optional
            `indexer_kwargs`, those slices will be applied to labels as well as
            posterior distribution parameters
        use_stats: Union[str, Collection[str]]
            Names of test statistics to perform.
            Must be in `constants.STATS_AVAILABLE`.
        drop: bool
            If True, drop quantiles when no quantiles requested or when 0
            psisamples taken.
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case testing is done
            over the remaining prefixes.
            All input dimensions are expected in distribution parameters.
            Slice will be applied to `labels` whenever the dimension is found
            in `labels`.

        Returns
        -------
        xr.Dataset
            With variable `pvalue` for statistics on approximate posterior means
            and `pvalue_quantiles` for quantiles of statistics on psisamples,
            if not dropped
        """
        a = self._approximate_alpha_core_prefix
        b = self._approximate_beta_core_prefix
        if indexer_kwargs:
            a = a.isel(indexer_kwargs)
            b = b.isel(indexer_kwargs)
            labels = labels.isel(indexer_kwargs, missing_dims="ignore")
        return _compute_stats(
            a,
            b,
            labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=drop,
        )


class MixinBootstrapPsi(MixinApproximatePsi, ABC):
    """Methods for PSI inference when bootstrap_alpha, bootstrap_beta defined"""

    @property
    @abstractmethod
    def bootstrap_alpha(self) -> xr.DataArray:
        """array(..., bootstrap_replicate) of alpha for posterior bootstrap replicates"""
        ...

    @property
    @abstractmethod
    def bootstrap_beta(self) -> xr.DataArray:
        """array(..., bootstrap_replicate) of beta for posterior bootstrap replicates"""
        ...

    @cached_property
    def _approximate_params(self) -> Tuple[xr.DataArray, xr.DataArray]:
        """Beta distribution parameters matching mean, variance of bootstrap mixture

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        alpha, beta = xr.apply_ufunc(
            bm.approximation,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[[], []],
            dask="allowed",
        )
        return (alpha, beta)

    @property
    def approximate_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx) alpha parameter of approximated bootstrap posterior

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        return self._approximate_params[0]

    @property
    def approximate_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx) beta parameter of approximated bootstrap posterior

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        return self._approximate_params[1]

    @cached_property
    def _bootstrap_moments(self) -> Tuple[xr.DataArray, xr.DataArray]:
        """Compute mean, variance of bootstrap posterior mixture"""
        agg_mean, agg_variance = xr.apply_ufunc(
            bm.moments,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[[], []],
            dask="allowed",
        )
        return (agg_mean, agg_variance)

    @property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        """array(...) means of mixtures of bootstrapped posteriors

        array(...) means of mixture of bootstrapped posterior distribution on PSI
        """
        return self._bootstrap_moments[0]

    @property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of mixtures of bootstrapped posteriors

        array(...) variances of mixtures of bootstrapped posterior
        distributions on PSI
        """
        return self._bootstrap_moments[1]

    @cached_property
    def bootstrap_psi_mean_legacy(self) -> xr.DataArray:
        """array(...) median of means of bootstrapped posteriors

        array(...) median of means of bootstrapped posterior distributions on PSI.

        Notes
        -----
        This is what was reported in MAJIQ v1 and v2.
        We have observed that if we increase the number of bootstrap replicates,
        both estimates tend close (but not exactly) to the raw posterior mean,
        which we now prefer.
        """
        return xr.apply_ufunc(
            bm.means_median,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            dask="allowed",
        )

    def bootstrap_cdf(
        self,
        x: Union[xr.DataArray, Sequence[float]],
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute cdf of mixture of bootstrapped posterior distribution

        Parameters
        ----------
        x: Union[xr.DataArray, Sequence[float]]
            potential realizations of distribution to evaluate CDF at
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that you *can* slice bootstrap_replicate, in which case,
            computation will be performed on the selected replicates.

        Returns
        -------
        xr.DataArray
            Return array(...) of CDF probabilities per connection. If
            `x` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "x"
        """
        alpha = self.bootstrap_alpha
        beta = self.bootstrap_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_cdf(alpha, beta, x)

    def bootstrap_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = (0.1, 0.9),
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute quantiles of mixture of bootstrapped posterior distributions

        Parameters
        ----------
        quantiles: Union[xr.DataArray, Sequence[float]]
            quantiles of distribution to compute
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that you *can* slice bootstrap_replicate, in which case,
            computation will be performed on the selected replicates.

        Returns
        -------
        xr.DataArray
            Return array(...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"

        Notes
        -----
        Please use `approximate_quantile` instead, which is faster, and what we
        think is a better representation of PSI
        """
        alpha = self.bootstrap_alpha
        beta = self.bootstrap_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_quantile(alpha, beta, quantiles=quantiles)

    def bootstrap_discretized_pmf(
        self,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        midpoint_approximation: bool = True,
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """Compute discretized PMF of bootstrap posterior mixture

        Parameters
        ----------
        nbins: int
            Number of uniform bins on [0, 1] on which probability mass will be
            computed
        midpoint_approximation: bool
            If midpoint_approximation, approximate PMF by computing
            unnormalized PDF over midpoints of uniform bins, then normalizing.
            Otherwise, use incomplete beta function.
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that you *can* slice bootstrap_replicate, in which case,
            computation will be performed on the selected replicates.
        """
        alpha = self.bootstrap_alpha
        beta = self.bootstrap_beta
        if indexer_kwargs:
            alpha = alpha.isel(indexer_kwargs)
            beta = beta.isel(indexer_kwargs)
        return _compute_posterior_discretized_pmf(
            alpha, beta, midpoint_approximation=midpoint_approximation, nbins=nbins
        )

    @cached_property
    def _bootstrap_alpha_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.bootstrap_alpha
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def _bootstrap_beta_core_prefix(self) -> xr.DataArray:
        """For computing statistics over a population of samples"""
        result = self.bootstrap_beta
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    def bootstrap_stats(
        self,
        labels: xr.DataArray,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        drop: bool = True,
        **indexer_kwargs: slice,
    ) -> xr.Dataset:
        """Statistics on bootstrap posterior means and psisamples.

        Statistics on bootstrap posterior means and quantiles of statistics
        from repeated sampling (psisamples) with respect to `labels`.

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles of pvalues to extract from repeating test statistics on
            `psisamples` draws from the posterior distributions
        psisamples: int
            Number of repeated draws from the posterior distributions to
            perform test statistics on
        labels: xr.DataArray
            Boolean labels over prefixes to perform tests on.
            Note that if it shares any dimensions with optional
            `indexer_kwargs`, those slices will be applied to labels as well as
            posterior distribution parameters
        use_stats: Union[str, Collection[str]]
            Names of test statistics to perform.
            Must be in `constants.STATS_AVAILABLE`.
        drop: bool
            If True, drop quantiles when no quantiles requested or when 0
            psisamples taken.
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case testing is done
            over the remaining prefixes.
            All input dimensions are expected in distribution parameters.
            Slice will be applied to `labels` whenever the dimension is found
            in `labels`.

        Returns
        -------
        xr.Dataset
            With variable `pvalue` for statistics on bootstrap posterior means
            and `pvalue_quantiles` for quantiles of statistics on psisamples,
            if not dropped
        """
        a = self._bootstrap_alpha_core_prefix
        b = self._bootstrap_beta_core_prefix
        if indexer_kwargs:
            a = a.isel(indexer_kwargs)
            b = b.isel(indexer_kwargs)
            labels = labels.isel(indexer_kwargs, missing_dims="ignore")
        return _compute_stats(
            a,
            b,
            labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=drop,
        )


class MixinRawTotalPopulation(ABC):
    """Methods for summarizing raw_total over a population (dimension prefix)"""

    @property
    @abstractmethod
    def raw_total(self):
        """array(...) total of coverage for event"""
        ...

    @property
    @abstractmethod
    def event_passed(self) -> xr.DataArray:
        """array(ec_idx, prefix) indicating if event passed for prefix"""
        ...

    @cached_property
    def _raw_total_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.raw_total.where(self.event_passed)
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def raw_total_population_median(self) -> xr.DataArray:
        """array(...) median over prefixes of `raw_total`"""
        return xr.apply_ufunc(
            nanmedian,
            self._raw_total_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    def raw_total_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of `raw_total`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case the quantile is
            taken over the remaining prefixes.

        Returns
        -------
        xr.DataArray
            array(..., `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        total = self._raw_total_core_prefix
        if indexer_kwargs:
            total = total.isel(indexer_kwargs)
        return _compute_population_quantile(
            total,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )


class MixinRawPsiMeanPopulation(ABC):
    """Methods for summarizing raw_psi_mean over a population (dimension prefix)"""

    @property
    @abstractmethod
    def raw_psi_mean(self):
        """array(...) means of raw posterior distribution on PSI"""
        ...

    @cached_property
    def _raw_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.raw_psi_mean
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def raw_psi_mean_population_median(self) -> xr.DataArray:
        """array(...) median over prefixes of `raw_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._raw_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    def raw_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of `raw_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case the quantile is
            taken over the remaining prefixes.

        Returns
        -------
        xr.DataArray
            array(..., `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        psi = self._raw_psi_mean_core_prefix
        if indexer_kwargs:
            psi = psi.isel(indexer_kwargs)
        return _compute_population_quantile(
            psi,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )


class MixinBootstrapPsiMeanPopulation(ABC):
    """Methods for summarizing bootstrap_psi_mean over a population (dimension prefix)"""

    @property
    @abstractmethod
    def bootstrap_psi_mean(self):
        """array(...) means of bootstrap posterior distribution on PSI"""
        ...

    @cached_property
    def _bootstrap_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.bootstrap_psi_mean
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": -1})
        return result

    @cached_property
    def bootstrap_psi_mean_population_median(self) -> xr.DataArray:
        """array(...) median over prefixes of `bootstrap_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._bootstrap_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    def bootstrap_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
        **indexer_kwargs: slice,
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of `bootstrap_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`
        indexer_kwargs: {dim: slice}, optional
            pairs of dimension names slices passed into `isel()` method on
            underlying xarray objects.
            Note that this can slice on "prefix", in which case the quantile is
            taken over the remaining prefixes.

        Returns
        -------
        xr.DataArray
            array(..., `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        psi = self._bootstrap_psi_mean_core_prefix
        if indexer_kwargs:
            psi = psi.isel(indexer_kwargs)
        return _compute_population_quantile(
            psi,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )


SelfT = TypeVar("SelfT", bound="MixinPsiOverPrefixes")


class MixinPsiOverPrefixes(
    MixinSubsettablePrefixes,
    MixinHasEvents,
    MixinRawPsi,
    MixinApproximatePsi,
    MixinRawPsiMeanPopulation,
    MixinBootstrapPsiMeanPopulation,
    MixinRawTotalPopulation,
):
    DIMS_BEFORE_PREFIX = ("ec_idx",)

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Initialize class with information on PSI over prefixes (`df`) and `events`"""
        # save df
        MixinSubsettablePrefixes.__init__(self, df)
        # save events
        MixinHasEvents.__init__(self, df, events)
        return

    @property
    def _non_df_args(self) -> Tuple:
        """Subsetting over prefixes requires positional argument about events"""
        return (self.events_df,)

    @property
    def _repr_other_dims_info(self) -> str:
        """Information about non-prefix dimensions for __repr__"""
        return f"{self.num_connections}"

    @abstractmethod
    def mask_events(self: SelfT, passed: xr.DataArray) -> SelfT:
        """Return class passing only events passed in input

        Parameters
        ----------
        passed: xr.DataArray
            boolean array(ec_idx) where connections marked as not passed
            (False) will be considered not passed over all prefixes
        """
        ...

    @property
    @abstractmethod
    def event_passed(self) -> xr.DataArray:
        """array(ec_idx, prefix) indicating if event passed for prefix"""
        ...

    @cached_property
    def num_passed(self) -> xr.DataArray:
        """Number of prefixes for which an event was passed"""
        return self.event_passed.sum("prefix")

    def passed_min_experiments(
        self,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> xr.DataArray:
        """Return boolean mask array for events passing min_experiments

        Parameters
        ----------
        min_experiments_f: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of prefixes that must pass
            individually for the event to be considered as passed for the group

        Returns
        -------
        xr.DataArray
            array: indicate whether an event passed min_experiments
        """
        return self.num_passed >= min_experiments(min_experiments_f, self.num_prefixes)

    def plot_violins(
        self,
        ec_idx: int,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        midpoint_approximation: bool = True,
        prefix_groups: Optional[Union[str, Mapping[str, str]]] = None,
        group_order: Optional[Sequence[str]] = None,
        ax=None,
        cmap=None,
        use_width: float = 0.9,
        rotate_group_labels: int = 45,
        swarm: bool = True,
        color_idx: Optional[int] = None,
    ) -> None:
        """Plot posterior distributions over groups of prefixes

        Parameters
        ----------
        ec_idx: int
            index of event connection to plot violins for
        nbins: int
            Compute PDF over endpoints of uniformly spaced bins on [0, 1].
            (the first and last values are computed at the midpoints in order
            to handle singularities at {0, 1} when either of the beta
            distribution parameters are less than 1).
        prefix_groups: Optional[Union[str, Mapping[str, str]]]
            If None, plot each prefix separately.
            If string, plot each prefix as a single group, with group_prefixes
            as its name.
            Otherwise, assumed to be a mapping from prefixes to group names
            (e.g. Dict, Series, etc.).
        group_order: Optional[Sequence[str]]
            If specified, the groups in the order they should be plotted.
            Otherwise, use the available groups in alphabetical order.
        use_width: float
            How much of plotting width each violin takes
        ax: Optional[Axes]
            Matplotlib axes to build plot in
        cmap: colors, optional
            matplotlib colors. If not provided, will use tab10 from matplotlib
        midpoint_approximation: bool
            Calculate PMF using approximation with midpoints of bins rather
            than incomplete beta function
        swarm: bool
            If True, use sns.swarmplot to plot PSI posterior means.
            If False, use sns.stripplot to plot PSI posterior means.
        color_idx: Optional[int]
            If None, use colors per group plotted.
            Otherwise, use specified index of palette as color.
        """
        import seaborn as sns

        if ax is None:
            try:
                import matplotlib.pyplot as plt
            except ImportError as err:
                raise ValueError(
                    "approximate_posterior_plot requires matplotlib"
                ) from err
            ax = plt.gca()
        if not (0 < use_width < 1):
            raise ValueError("use_width must be between 0 and 1")
        if cmap is None:
            cmap = sns.color_palette()  # default seaborn color palette

        def get_color(idx):
            try:
                return cmap[idx]
            except TypeError:
                return cmap(idx)

        pmf = (
            self.approximate_discretized_pmf(
                nbins=nbins,
                midpoint_approximation=midpoint_approximation,
                ec_idx=slice(ec_idx, 1 + ec_idx),
            )
            .isel(ec_idx=0)
            .load()
        )
        df_psi = self.raw_psi_mean.isel(ec_idx=ec_idx).to_dataframe("PSI")

        if prefix_groups is None:
            df_psi = df_psi.assign(group=df_psi.index)
        elif isinstance(prefix_groups, str):
            df_psi = df_psi.assign(group=prefix_groups)
        else:
            df_psi = df_psi.assign(group=[prefix_groups[x] for x in df_psi.index])
        dfg_psi = df_psi.groupby("group")
        group_order = group_order or sorted(dfg_psi.groups)
        # override colors if color_idx
        if color_idx is None:
            color_kwargs = {"palette": cmap}
        else:
            color_kwargs = {"color": get_color(color_idx)}
        psi_mean_plot = sns.swarmplot if swarm else sns.stripplot
        psi_mean_plot(
            data=df_psi, x="group", y="PSI", order=group_order, ax=ax, **color_kwargs
        )
        # plot violins on top of the swarmplots
        for offset, group in enumerate(group_order):
            # prefixes for this group
            group_prefixes = dfg_psi.groups[group]
            # color for this group
            color = get_color(offset if color_idx is None else color_idx)
            # pmf over group
            pmf_group = (
                # get PMF for prefixes in group, averaging out everything but bins
                pmf.sel(prefix=group_prefixes)
                .mean(list(set(pmf.dims) - {"pmf_bin"}))
                # normalize to width for plotting
                .pipe(lambda x: x * 0.5 * use_width / x.max())
                # assign dilated coordinates from start/end bins
                .assign_coords(plot_psi=("pmf_bin", np.linspace(0, 1, nbins)))
                # ignore negligible values of p
                .pipe(lambda x: x[x > 5e-4 * use_width])
            )
            ax.fill_betweenx(
                pmf_group.plot_psi,
                offset - pmf_group,
                offset + pmf_group,
                fc=color,
                ec=color,
                alpha=0.7,
            )
        ax.set_xticklabels(
            ax.get_xticklabels(), rotation=rotate_group_labels, va="top", ha="right"
        )
        ax.set_ylim(0, 1)
        return


def _compute_population_quantile(
    x: xr.DataArray,
    quantiles: Sequence[float],
    quantile_dim_name: str = "population_quantile",
) -> xr.DataArray:
    """Compute quantiles (ignoring nan) over dimension prefix"""
    quantiles_xr = xr.DataArray(quantiles, [(quantile_dim_name, quantiles)])
    return xr.apply_ufunc(
        nanquantile,
        x,
        quantiles_xr,
        input_core_dims=[["prefix"], [quantile_dim_name]],
        output_core_dims=[[quantile_dim_name]],
        dask="allowed",
    )


def _compute_stats(
    a: xr.DataArray,
    b: xr.DataArray,
    labels: xr.DataArray,
    quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
    psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
    use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    mix_dim: str = "bootstrap_replicate",
    drop: bool = True,
) -> xr.Dataset:
    """Get pvalues on distribution means and pvalue quantiles on psisamples"""
    # make sure they have mixture dimension
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    # define stats vector to input
    if isinstance(use_stats, str):
        use_stats = [use_stats]  # make sure it's always a collection
    else:
        use_stats = sorted(set(use_stats))
    if any(unknown := [x for x in use_stats if x not in constants.STATS_AVAILABLE]):
        raise ValueError(
            f"Input statistics {unknown} unavailable"
            f" (available {set(constants.STATS_AVAILABLE.keys())})"
        )
    use_stats_xr: xr.DataArray = xr.DataArray(
        [constants.STATS_AVAILABLE[x] for x in use_stats], [("stats", use_stats)]
    )
    # quantiles to input
    quantiles_xr = xr.DataArray(
        _ := np.array(quantiles, dtype=np.float64), [("pval_quantile", _)]
    )
    # get result from stats_sample
    result = xr.apply_ufunc(
        bm.stats_sample,
        a,
        b,
        labels,
        quantiles_xr,
        psisamples,
        use_stats_xr,
        input_core_dims=[
            ["prefix", mix_dim],
            ["prefix", mix_dim],
            ["prefix"],
            ["pval_quantile"],
            [],
            ["stats"],
        ],
        output_core_dims=[["stats"], ["stats", "pval_quantile"]],
        output_dtypes=[np.float64, np.float64],
        dask="allowed",
    )
    # convert result into xr.Dataset
    ds = xr.Dataset(
        {
            "pvalue": result[0],
            "pvalue_quantiles": result[1].assign_attrs(psisamples=psisamples),
        },
    )
    # drop psisample quantiles if none taken and drop requested
    if drop and not (psisamples > 0 and quantiles):
        ds = ds.drop_dims("pval_quantile")
    return ds


def _compute_posterior_cdf(
    a: xr.DataArray,
    b: xr.DataArray,
    x: Union[xr.DataArray, Sequence[float]],
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Caclulate cdf over posterior distribution

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
    x: Union[xr.DataArray, Sequence[float]]
        potential realizations of the distribution at which to evaluate the CDF
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution
    """
    if not isinstance(x, xr.DataArray):
        x_arr = np.array(x, dtype=a.dtype)
        if x_arr.ndim > 1:
            raise ValueError("Unable to handle non-xarray multi-dimensional x")
        elif x_arr.ndim == 0:
            x_arr = x_arr[np.newaxis]
        x = xr.DataArray(x_arr, [("x", x_arr)])
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    if a.sizes[mix_dim] != b.sizes[mix_dim]:
        raise ValueError(
            f"a and b must share same size for dimension {mix_dim}"
            f" ({a.sizes = }, {b.sizes = })"
        )
    return xr.apply_ufunc(
        bm.cdf,
        x,
        a,
        b,
        input_core_dims=[[], [mix_dim], [mix_dim]],
        dask="allowed",
    )


def _compute_posterior_quantile(
    a: xr.DataArray,
    b: xr.DataArray,
    quantiles: Union[xr.DataArray, Sequence[float]] = (0.1, 0.9),
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Calculate quantiles over posterior distribution

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
    quantiles: Union[xr.DataArray, Sequence[float]]
        quantiles to compute over posterior distribution. If Sequence[float],
        add "quantiles" dimension to result for input quantiles
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution
    """
    if not isinstance(quantiles, xr.DataArray):
        quantiles_arr = np.array(quantiles, dtype=a.dtype)
        if quantiles_arr.ndim > 1:
            raise ValueError("Unable to handle non-xarray multi-dimensional quantiles")
        elif quantiles_arr.ndim == 0:
            quantiles_arr = quantiles_arr[np.newaxis]
        quantiles = xr.DataArray(quantiles_arr, [("quantiles", quantiles_arr)])
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    if a.sizes[mix_dim] != b.sizes[mix_dim]:
        raise ValueError(
            f"a and b must share same size for dimension {mix_dim}"
            f" ({a.sizes = }, {b.sizes = })"
        )
    if quantiles.size and a.size and b.size:
        return xr.apply_ufunc(
            bm.quantile,
            quantiles,
            a,
            b,
            input_core_dims=[[], [mix_dim], [mix_dim]],
            dask="allowed",
        )
    else:
        # if the result will have size zero, bm.quantile sometimes raises
        # ValueError: operands could not be broadcast together with remapped
        # shapes. We can get desired behavior by broadcasting/reducing using
        # xr.dot
        return xr.dot(quantiles, a, b, dim=[mix_dim])


def _get_pmf_bins(nbins: int, dtype=np.float32) -> xr.DataArray:
    """Array over dimension "pmf_bin" for `nbins` bins over [0, 1]"""
    endpoints = np.linspace(0, 1, 1 + nbins, dtype=dtype)
    midpoints = endpoints[:-1] + 0.5 * np.diff(endpoints)
    return xr.DataArray(
        midpoints,
        {
            "pmf_bin_start": ("pmf_bin", endpoints[:-1]),
            "pmf_bin_end": ("pmf_bin", endpoints[1:]),
            "pmf_bin_mid": ("pmf_bin", midpoints),
        },
        dims=["pmf_bin"],
    )


def _compute_posterior_discretized_pmf(
    a: xr.DataArray,
    b: xr.DataArray,
    midpoint_approximation: bool = True,
    nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    if midpoint_approximation:
        return _compute_posterior_discretized_pmf_approximation(
            a, b, nbins=nbins, mix_dim=mix_dim
        )
    else:
        return _compute_posterior_discretized_pmf_exact(
            a, b, nbins=nbins, mix_dim=mix_dim
        )


def _compute_posterior_discretized_pmf_exact(
    a: xr.DataArray,
    b: xr.DataArray,
    nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Compute discretized PMF of posterior (mixture)

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
        Dimensions: [..., ?mix_dim]
    nbins: int
        Number of uniform bins on [0, 1] on which probability mass will be
        computed
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution

    Returns
    -------
    xr.DataArray
        discretized PMF for posterior (dimensions: [..., pmf_bin])
    """
    dummy_bins = _get_pmf_bins(nbins, dtype=a.dtype)
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    return xr.apply_ufunc(
        bm.pmf,
        a,
        b,
        dummy_bins,
        input_core_dims=[
            [mix_dim],
            [mix_dim],
            ["pmf_bin"],
        ],
        output_core_dims=[["pmf_bin"]],
        dask="allowed",
    )


def _compute_posterior_discretized_pmf_approximation(
    a: xr.DataArray,
    b: xr.DataArray,
    nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Approximation of discretized PMF of posterior (mixture)

    Approximation of discretized PMF of posterior (mixture).
    Computes unnormalized PDF at midpoints, then normalizes over bins.
    If there is a mixture dimension, takes the average over that dimension.

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
        Dimensions: [..., ?mix_dim]
    nbins: int
        Number of uniform bins on [0, 1] on which probability mass will be
        computed
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution

    Returns
    -------
    xr.DataArray
        discretized PMF for posterior (dimensions: [..., pmf_bin])

    Notes
    -----
    This is appropriate for plotting because it is much faster than computing
    the PMF and is close to the true PMF.

    This is closely related to using a midpoint Riemann sum to compute a Beta
    function.
    """
    # points to evaluate at
    x = _get_pmf_bins(nbins, dtype=a.dtype)
    # log of the unnormalized PDF is
    log_updf = (a - 1) * np.log(x) + (b - 1) * np.log1p(-x)
    # log normalization constant over bins
    log_updf_Z = xr.apply_ufunc(
        _logsumexp, log_updf, input_core_dims=[["pmf_bin"]], dask="allowed"
    )
    # do we have the mixture dimension?
    if mix_dim not in log_updf.dims:
        return np.exp(log_updf - log_updf_Z)
    # otherwise, we need to take the average over the mixture dimension
    log_pmf_sum = xr.apply_ufunc(
        _logsumexp, log_updf - log_updf_Z, input_core_dims=[[mix_dim]], dask="allowed"
    )
    # turn to proper probability and normalize over bins
    pmf_sum = np.exp(log_pmf_sum)
    return pmf_sum / pmf_sum.sum("pmf_bin")


def _logsumexp(x):
    """Do logsumexp over the last axis of x

    Adapted from scipy.special.logsumexp, which is not used because it eagerly
    computes on dask arrays.
    """
    x_max = np.max(x, axis=-1)
    # TODO handle nan values, infinite values?
    exp_x_scaled = np.exp(x - x_max[..., np.newaxis])
    sumexp_x_scaled = np.sum(exp_x_scaled, axis=-1)
    logsumexp_x_scaled = np.log(sumexp_x_scaled)
    return logsumexp_x_scaled + x_max
