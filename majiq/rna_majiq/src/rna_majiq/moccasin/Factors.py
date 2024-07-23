"""
Factors.py

Observed and/or modeled confounding or non-confounding factors for modeling
PsiCoverage

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Collection, Dict, Final, List, Optional, Tuple, Union

import numpy as np
import xarray as xr


class Factors(object):
    """Numeric factors over prefixes for modeling :class:`PsiCoverage`

    Parameters
    ----------
    df: xr.Dataset
        With variables as in `Factors.EXPECTED_VARIABLES`.
    """

    EXPECTED_VARIABLES: Final[Dict[str, Tuple[str, ...]]] = {
        "prefix": ("prefix",),
        "factor": ("factor",),
        "confounding": ("factor",),
        "factors": ("prefix", "factor"),
    }
    DTYPE = np.float32
    # ZARR_GROUP = "moccasin_factors"
    ZARR_GROUP = None

    def __init__(self, df: xr.Dataset):
        # save df after checking validity
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        for v in df.variables.values():
            v.encoding.clear()
        self.df: Final[xr.Dataset] = df.load()
        if self.df["factors"].isnull().any():
            raise ValueError("Cannot have nan-valued factors")
        return

    def __repr__(self) -> str:
        return (
            f"Factors[{self.num_prefixes} prefixes, {self.num_factors} factors]"
            f" ({self.num_confounding} confounding)"
        )

    def factor_report(self) -> str:
        """Two column table enumerating name of each factor and if confounding"""
        return "\n".join(
            f"Factor {factor_name}\t{factor_type}"
            for factor_name, factor_type in self.confounding.to_series()
            .map({False: "nonconfounding", True: "confounding"})
            .items()
        )

    @property
    def factors(self) -> xr.DataArray:
        """array with factors for each prefix"""
        return self.df["factors"]

    @property
    def num_factors(self) -> int:
        return self.df.sizes["factor"]

    @property
    def confounding(self) -> xr.DataArray:
        """array indicating if a factor is confounding or not"""
        return self.df["confounding"]

    @property
    def num_confounding(self) -> int:
        return self.confounding.sum().values[()]

    @property
    def factor_names(self) -> List[str]:
        return self.df["factor"].values.tolist()

    @property
    def prefixes(self) -> List[str]:
        return self.df["prefix"].values.tolist()

    @property
    def num_prefixes(self) -> int:
        return self.df.sizes["prefix"]

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Save :class:`Factors` to specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        """
        save_df = self.df.chunk(self.df.sizes)
        save_df.to_zarr(
            path,
            mode="w",
            group=self.ZARR_GROUP,
        )
        return

    def combine(self, add_factors: "Factors") -> "Factors":
        """Create updated :class:`Factors` with `new_factors`

        Requires either dimension "factor" or "new_factor" in `new_factors`
        """
        if unique := set(self.prefixes).symmetric_difference(set(add_factors.prefixes)):
            raise ValueError(
                f"Combined factors require identical prefixes ({unique = })"
            )
        if overlap := set(self.factor_names) & set(add_factors.factor_names):
            raise ValueError(
                f"Combined factors require unique factor names ({overlap = })"
            )
        return Factors(xr.concat([self.df, add_factors.df], dim="factor"))

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path, List[Union[str, Path]]],
        prefixes: Optional[Union[str, List[str]]] = None,
    ) -> "Factors":
        """Load :class:`Factors` from one or more specified Zarr paths"""
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=cls.ZARR_GROUP,
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="override",  # use first file if overlapping prefixes
            coords="minimal",
            data_vars="minimal",
        ).drop_duplicates("prefix")
        # load factors for desired prefixes into memory
        if prefixes:
            if isinstance(prefixes, str):
                prefixes = [prefixes]
            df = df.sel(prefix=prefixes)
        df = df
        return Factors(df)

    @classmethod
    def from_tsv(
        cls,
        factors_tsv: Union[str, Path, List[Union[str, Path]]],
        confounding: Union[str, Collection[str]],
        prefixes: Optional[Union[str, List[str]]] = None,
    ) -> "Factors":
        """Load :class:`Factors` from `factors_tsv`, marking `confounding` variables

        Parameters
        ----------
        factors_tsv: Union[str, Path]
            Path for TSV with model matrix for all prefixes being processed.
            Required column: "prefix".
        confounding: Union[str, Collection[str]]
            Names of confounding variables in `factors_tsv`. All others will be
            considered nonconfounding.
        prefixes: Optional[Union[str, List[str]]]
            If specified, subset of prefixes that will specifically used.
        """
        import pandas as pd

        if not isinstance(factors_tsv, list):
            factors_tsv = [factors_tsv]

        dfs = [
            pd.read_csv(x, sep="\t")
            .astype({"prefix": str})
            .set_index("prefix", verify_integrity=True)
            for x in factors_tsv
        ]
        if not all(df.columns.symmetric_difference(dfs[0].columns).empty for df in dfs):
            raise ValueError(
                "Input TSV files do not have exact same columns\n"
                + "\n".join(
                    f"{path}: {df.columns}" for path, df in zip(factors_tsv, dfs)
                )
            )

        df = pd.concat(dfs, join="inner")
        # drop duplicate prefixes, if any
        df = df.loc[~df.index.duplicated(keep="first")]

        if prefixes:
            if isinstance(prefixes, str):
                prefixes = [prefixes]
            df = df.loc[prefixes]
        df = df.astype(cls.DTYPE)
        if isinstance(confounding, str):
            confounding = [confounding]
        if missing := set(x for x in confounding if x not in df.columns):
            raise ValueError(
                f"Not all specified confounders were found from TSV ({missing = })"
            )
        ds = xr.Dataset(
            {
                "factors": (("prefix", "factor"), df.values),
            },
            {
                "prefix": df.index,
                "factor": df.columns,
                "confounding": ("factor", [x in confounding for x in df.columns]),
            },
        )
        return Factors(ds)

    @classmethod
    def intercept_only(cls, prefixes: Union[str, List[str]]) -> "Factors":
        """Create :class:`Factors` for `prefixes` with nonconfounding intercept term"""
        df = xr.Dataset(
            {
                "factors": (
                    ("prefix", "factor"),
                    np.ones((len(prefixes), 1), dtype=cls.DTYPE),
                ),
            },
            {
                "prefix": prefixes,
                "factor": ["intercept"],
                "confounding": ("factor", [False]),
            },
        )
        return Factors(df)

    def __getitem__(self, prefixes) -> "Factors":
        """Subset :class:`Factors` to selected prefixes"""
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        return Factors(self.df.sel(prefix=prefixes))
