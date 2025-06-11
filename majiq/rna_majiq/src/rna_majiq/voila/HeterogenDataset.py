"""
HeterogenDataset.py

Intermediate quantifications for HET for use in VOILA visualizations

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Collection, Dict, Final, List, Optional, Sequence, Union, cast

import dask
import numpy as np
import pandas as pd
import xarray as xr
from dask.delayed import Delayed
from moccasin.moccasin import persist_with_progress

import rna_majiq.constants as constants

from ..core.Heterogen import Heterogen
from ..core.MixinHasEvents import MixinHasEvents
from ..core.PsiGroup import PsiGroup
from ..core.SpliceGraph import SpliceGraph


class HeterogenDataset(MixinHasEvents):
    """Precomputed Heterogen quantifications for use in VOILA visualizations"""

    EXPECTED_VARIABLES: Final = {
        # identity of comparisons
        "comparison_grp1": ("comparison",),
        "comparison_grp2": ("comparison",),
        # identity of groups
        "grp": ("grp",),
        "grp_size": ("grp",),  # check vs prefix_grp
        # identity of prefixes
        "prefix": ("prefix",),
        "prefix_grp": ("prefix",),
        # identity of stats
        "stats": ("stats",),
        # psi for prefixes
        "raw_alpha": ("ec_idx", "prefix"),
        "raw_beta": ("ec_idx", "prefix"),
        "approximate_alpha": ("ec_idx", "prefix"),
        "approximate_beta": ("ec_idx", "prefix"),
        # psisamples for comparisons
        "psisamples": ("comparison",),
        "pval_quantile": ("pval_quantile",),
        # stats
        "raw_pvalue": ("comparison", "ec_idx", "stats"),
        "approximate_pvalue": ("comparison", "ec_idx", "stats"),
        "approximate_pvalue_quantiles": (
            "comparison",
            "ec_idx",
            "stats",
            "pval_quantile",
        ),
    }

    @classmethod
    def from_heterogen(
        cls,
        het: Heterogen,
        pvalue_quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
    ) -> "HeterogenDataset":
        if prefix_overlap := set(het.psi1.prefixes) & set(het.psi2.prefixes):
            # do not allow HeterogenDataset with overlapping prefixes
            raise ValueError(
                "HeterogenDataset cannot have prefixes belonging to both groups"
                f" (overlapping prefixes = {prefix_overlap}"
            )
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
        # compute stats
        raw_stats = het.raw_stats(use_stats=use_stats).expand_dims(comparison=1)
        approximate_stats = het.approximate_stats(
            quantiles=pvalue_quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=False,
        ).expand_dims(comparison=1)
        # parameters for alpha/beta of raw/approximate distributions
        raw_alpha = xr.concat([het.psi1.raw_alpha, het.psi2.raw_alpha], dim="prefix")
        raw_beta = xr.concat([het.psi1.raw_beta, het.psi2.raw_beta], dim="prefix")
        approximate_alpha = xr.concat(
            [het.psi1.approximate_alpha, het.psi2.approximate_alpha], dim="prefix"
        )
        approximate_beta = xr.concat(
            [het.psi1.approximate_beta, het.psi2.approximate_beta], dim="prefix"
        )
        # construct dataset
        df = xr.Dataset(
            {
                "raw_alpha": raw_alpha,
                "raw_beta": raw_beta,
                "approximate_alpha": approximate_alpha,
                "approximate_beta": approximate_beta,
                "raw_pvalue": raw_stats["pvalue"],
                "approximate_pvalue": approximate_stats["pvalue"],
                "approximate_pvalue_quantiles": approximate_stats["pvalue_quantiles"],
            },
            {
                "comparison_grp1": ("comparison", [het.name1]),
                "comparison_grp2": ("comparison", [het.name2]),
                "grp": [het.name1, het.name2],
                "grp_size": ("grp", [het.psi1.num_prefixes, het.psi2.num_prefixes]),
                # "prefix": het.psi1.prefixes + het.psi2.prefixes,
                "prefix_grp": (
                    "prefix",
                    [het.name1] * het.psi1.num_prefixes
                    + [het.name2] * het.psi2.num_prefixes,
                ),
                # "stats": use_stats,
                "psisamples": ("comparison", [psisamples]),
            },
            {
                "prefixes_grp1": het.psi1_original_prefix,
                "prefixes_grp2": het.psi2_original_prefix,
            }
        )
        return HeterogenDataset(
            df,
            het.events_df,
            hold_temporary=[*het.psi1.hold_temporary, *het.psi2.hold_temporary],
        )

    def __init__(
        self,
        df: xr.Dataset,
        events: xr.Dataset,
        hold_temporary: Collection[TemporaryDirectory] = tuple(),
    ):
        """Initialize :class:`HeterogenDataset` with specified underlying data

        Parameters
        ----------
        df: xr.Dataset
            Variables/coordinates matching HeterogenDataset.EXPECTED_VARIABLES
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        hold_temporary: Optional[Collection[TemporaryDirectory]]
            Hold onto temporary directories, preventing their deletion (used by
            :meth:`PsiGroup.from_psicov` for temporary store of
            rechunked beta distribution parameters)
        """
        # save events
        MixinHasEvents.__init__(self, df, events)
        # save df appropriately after checking validity
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        # verify that all compared groups are recognized
        if not (
            df["comparison_grp1"].load().isin(df["grp"]).all()
            and df["comparison_grp2"].load().isin(df["grp"]).all()
        ):
            raise ValueError(
                "Not all compared groups"
                f" ({df['comparison_grp1'] = }, {df['comparison_grp2'] = })"
                f" have group information on dimension {df['grp'] = }"
            )
        # verify that prefixes/groups match up appropriately
        if (
            not df["prefix_grp"]
            .to_series()
            .value_counts()
            .sort_index()
            .equals(df["grp_size"].to_series().sort_index())
        ):
            raise ValueError("Mismatch between prefix_grp numbers and grp_size")
        # make these values available
        df["prefix_grp"].load()
        with dask.config.set(**{"array.slicing.split_large_chunks": False}):
            self.groups: Final[Dict[str, PsiGroup]] = {
                grp: PsiGroup(
                    df[
                        [
                            "raw_alpha",
                            "raw_beta",
                            "approximate_alpha",
                            "approximate_beta",
                        ]
                    ].sel(prefix=df["prefix_grp"] == grp),
                    self.events_df,
                )
                for grp in df["grp"].values
            }
        self.df: Final[xr.Dataset] = df.transpose(
            "comparison", "grp", "ec_idx", "prefix", "stats", "pval_quantile"
        ).chunk({"prefix": -1})
        self.hold_temporary: Final = hold_temporary
        return

    def to_zarr(
        self,
        path: Union[str, Path],
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save :class:`HeterogenDataset` to specified path"""
        save_df_future = cast(
            Delayed,
            self.df.to_zarr(
                path,
                mode="w",
                group=constants.NC_HETEROGEN,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = persist_with_progress(save_df_future)
        save_df_future.compute()
        self.events_to_zarr(path, mode="a", consolidated=consolidated)
        return

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path, List[Union[str, Path]]]
    ) -> "HeterogenDataset":
        """Load :class:`HeterogenDataset` from one or more specified paths"""
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_HETEROGEN,
            combine="nested",
            preprocess=lambda ds: ds.assign_coords(comparison=[ds.encoding["source"]]),
            join="outer",
            # no conflicts means that masked values (from not passing
            # comparison filters) from one file will be obtained from unmasked
            # values with same coordinates in another file (if it passed
            # comparison filters there)
            compat="no_conflicts",
            coords="minimal",
            data_vars="minimal",
        )
        # make sure grp_size is int (sometimes switched to float)
        # (probably because of compat="no_conflicts")
        df["grp_size"] = df["grp_size"].astype(np.int64)
        if len(path) > 2:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        return HeterogenDataset(df, events_df)

    @property
    def num_comparisons(self) -> int:
        """Number of comparisons in this dataset"""
        return self.df.sizes["comparison"]

    @property
    def num_groups(self) -> int:
        """Number of groups in this dataset"""
        return self.df.sizes["grp"]

    @property
    def num_prefixes(self) -> int:
        return self.df.sizes["prefix"]

    @property
    def comparisons(self):
        """Enumerate which groups were compared in this dataset"""
        return tuple(
            zip(self.df["comparison_grp1"].values, self.df["comparison_grp2"].values)
        )

    def __repr__(self) -> str:
        return (
            f"HeterogenDataset[{self.num_connections}]"
            f" for {self.num_comparisons} comparisons of"
            f" {self.num_groups} groups, {self.num_prefixes} total prefixes"
        )

    def psi_pmf(
        self,
        grp: str,
        ec_idx: Optional[slice],
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        midpoint_approximation: bool = True,
    ) -> xr.DataArray:
        """Average of PSI approximate bootstrap posterior distributions for group

        Parameters
        ----------
        grp: str
            quantification group to get average approximate distribution for
        ec_idx: Optional[slice]
            Slice of ec_idx to compute PMF for. If None explicitly passed,
            compute for all ec_idx (this may be an expensive computation)
        nbins: int
            Number of equally-sized bins on [0, 1] on which probability mass
            will be computed/averaged.
        midpoint_approximation: bool
            Set to false to get exact calculation rather than approximation
        """
        ec_select: Dict[str, slice] = dict()
        if ec_idx:
            ec_select["ec_idx"] = ec_idx
        return self.groups[grp].approximate_discretized_pmf(
            nbins=nbins, midpoint_approximation=midpoint_approximation, **ec_select
        )

    @property
    def psisamples(self) -> xr.DataArray:
        return self.df["psisamples"]

    @property
    def num_stats(self) -> int:
        return self.df.sizes["stats"]

    @property
    def stats(self) -> xr.DataArray:
        return self.df["stats"]

    @property
    def pval_quantile(self) -> xr.DataArray:
        return self.df["pval_quantile"]

    @property
    def raw_pvalue(self) -> xr.DataArray:
        return self.df["raw_pvalue"]

    @property
    def approximate_pvalue(self) -> xr.DataArray:
        return self.df["approximate_pvalue"]

    @property
    def approximate_pvalue_quantiles(self) -> xr.DataArray:
        return self.df["approximate_pvalue_quantiles"]

    @cached_property
    def any_passed(self) -> xr.DataArray:
        return self.df["raw_alpha"].notnull().any("prefix")

    @property
    def comparison_experiments(self) -> list:
        return [self.df.prefixes_grp1, self.df.prefixes_grp2]

    def to_dataframe(
        self,
        sg: Optional[SpliceGraph] = None,
        annotated: Optional[SpliceGraph] = None,
        population_quantiles: Sequence[
            float
        ] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        show_progress: bool = False,
    ) -> pd.DataFrame:
        """Return table of quantifications for TSV output

        Parameters
        ----------
        sg: Optional[SpliceGraph]
            If provided, splicegraph with introns/junctions consistent with
            events used to annotate resulting dataframe
        annotated: Optional[SpliceGraph]
            If specified, override denovo definitions by comparing
            exons/introns/junctions to annotated splicegraph.
            Only used if `sg` also provided.
        population_quantiles: Sequence[float]
            quantiles of raw_psi_mean to evaluate per group
        show_progress: bool
            show progress bar in dask if enabled
        """
        # build tables with columns that will be concatenated together
        concat_df: List[pd.DataFrame] = list()
        # population quantiles to compute (add median, unique up to 3 decimal places)
        population_quantiles = sorted({0.5, *np.round(population_quantiles, 3)})
        # add dataframe with events annotations
        if sg is not None:
            concat_df.append(
                self.get_events(sg.introns, sg.junctions).ec_dataframe(
                    annotated=annotated
                )
            )
        # build dataset with quantifications to load simultaneously
        ds: xr.Dataset = self.df.drop_dims("prefix").assign(
            any_passed=self.any_passed,
            raw_psi_quantile=xr.concat(
                [
                    x.raw_psi_mean_population_quantile(
                        population_quantiles
                    ).expand_dims(grp=[grp])
                    for grp, x in self.groups.items()
                ],
                dim="grp",
            ),
            num_passed=xr.concat(
                [x.num_passed.expand_dims(grp=[grp]) for grp, x in self.groups.items()],
                dim="grp",
            ),
        )
        # load/compute into memory
        if show_progress:
            ds = persist_with_progress(ds)
        ds = ds.load()

        # add number of experiments each event connection passed
        concat_df.append(
            # num_passed has dimensions ("grp", "ec_idx")
            ds["num_passed"]
            # so the series has a multiindex with grp, ec_idx
            .to_series()
            # unstack puts the levels from grp as columns
            .unstack("grp")
            # updates the column names to contextualize vs others in concat_df
            .add_suffix("-num_passed")
        )

        # add raw_psi_quantile into table
        df_quantile = (
            ds["raw_psi_quantile"]
            .to_series()
            .unstack(["population_quantile", "grp"])
            .sort_index(axis=1)
        )
        df_quantile.columns = [
            f"{grp}-raw_psi_quantile_{q:0.3f}" for q, grp in df_quantile.columns
        ]
        concat_df.append(df_quantile)
        del df_quantile

        # add comparisons/statistics into table
        for i, (grp1, grp2) in enumerate(self.comparisons):
            comparison_prefix = f"{grp1}-vs-{grp2}-" if self.num_comparisons > 1 else ""
            for j, statistic in enumerate(self.stats.values):
                stat_prefix = f"{statistic}-" if self.num_stats > 1 else ""
                concat_df.append(
                    ds[["raw_pvalue", "approximate_pvalue"]]
                    .isel(comparison=i, stats=j, drop=True)
                    .to_dataframe()
                    .assign(
                        **{
                            f"approximate_pvalue_quantiles_{q:0.3f}": ds[
                                "approximate_pvalue_quantiles"
                            ]
                            .isel(comparison=i, stats=j, pval_quantile=k, drop=True)
                            .to_series()
                            for k, q in enumerate(self.pval_quantile.values)
                        }
                    )
                    .add_prefix(f"{comparison_prefix}{stat_prefix}")
                )

        return pd.concat(concat_df, axis=1, join="inner").loc[ds["any_passed"].values]
