"""
MixinSubsettablePrefixes.py

Abstract mixin class to define objects with prefixes which hold onto df:
xr.Dataset with information on prefixes which can be subset to different
prefixes (e.g. PsiCoverage, PsiGroup, SpliceGraphReads).

Author: Joseph K Aicher
"""

from abc import ABC, abstractmethod
from typing import (
    ClassVar,
    Collection,
    Dict,
    Final,
    List,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
)

import numpy as np
import numpy.typing as npt
import xarray as xr

SelfT = TypeVar("SelfT", bound="MixinSubsettablePrefixes")


class MixinSubsettablePrefixes(ABC):
    """Abstract mixin class to define objects with prefixes

    Information over prefixes are stored in an xarray Dataset `df`.
    This should be subsettable over the "prefix" dimension.

    All deriving classes should have `df` as its first argument.
    In order to permit making subsets with derived classes that have more
    constructor arguments, derived classes must specialize the properties
    `_non_df_args` and `_non_df_kwargs` which will be passed into the
    constructor after `df`.

    Parameters
    ----------
    df: xr.Dataset
        xarray dataset with information over prefixes which can be subset.
        It must have variables/coordinates as defined in `EXPECTED_VARIABLES`
        (keys are names, values are dimensions they must have).
        The arrays will be transposed to have `DIMS_BEFORE_PREFIX` before the
        prefix dimension, with all other dimensions following afterwards.
    *_non_df_args, **_non_df_kwargs
        unused arguments for constructor that are used for derived classes
    """

    EXPECTED_VARIABLES: ClassVar[Dict[str, Tuple[str, ...]]] = {"prefix": ("prefix",)}
    DIMS_BEFORE_PREFIX: ClassVar[Sequence[str]] = tuple()

    def __init__(self, df: xr.Dataset, *_non_df_args, **_non_df_kwargs):
        """Initialize class with information over prefixes stored in df"""
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        # remove encoding
        for v in df.variables.values():
            v.encoding.clear()
        # transpose the data as appropriate
        self.df: Final[xr.Dataset] = df.transpose(
            *self.DIMS_BEFORE_PREFIX, "prefix", ...
        )
        return

    @property
    def _non_df_args(self) -> Tuple:
        """Additional positional arguments for constructor after `df`"""
        return tuple()

    @property
    def _non_df_kwargs(self) -> Dict:
        """Additional keyword arguments for constructor after `df`, `_non_df_args`"""
        return dict()

    def _with_updated_df(
        self: SelfT,
        df: xr.Dataset,
        override_args: Optional[Sequence] = None,
        update_kwargs: Optional[Dict] = None,
    ) -> SelfT:
        """Create updated instance of same type as `self` with updated `df`

        Parameters
        ----------
        df: xr.Dataset
            Use this xarray dataset for the class (first argument to
            constructor)
        override_args: Optional[Sequence]
            Instead of using `self._non_df_args` for additional positional
            arguments, use these positional arguments
        update_kwargs: Dict
            Add or replace elements of `self._non_df_kwargs` with values in
            input dictionary for keyword arguments to constructor
        """
        if update_kwargs is None:
            update_kwargs = {}
        extra_args = self._non_df_args if override_args is None else override_args
        extra_kwargs = {**self._non_df_kwargs, **update_kwargs}
        return self.__class__(df, *extra_args, **extra_kwargs)

    @property
    def num_prefixes(self) -> int:
        """Number of independent experiments"""
        return self.df.sizes["prefix"]

    @property
    def prefixes(self) -> List[str]:
        """Prefixes: Names of independent units of analysis

        Names of independent units of analysis (e.g. experiments, aggregations
        over experiments). Name derived from prefixes of BAM files used as
        experiment names.
        """
        return self.df["prefix"].values.tolist()

    @property
    @abstractmethod
    def _repr_other_dims_info(self) -> str:
        """Information about non-prefix dimensions for __repr__"""
        ...

    def __repr__(self) -> str:
        """String representation"""
        MAX_PREFIXES_END = 1  # how many prefixes on either end to display
        if self.num_prefixes > 2 * MAX_PREFIXES_END:
            print_prefixes_list = [
                *self.prefixes[:MAX_PREFIXES_END],
                *([] if self.num_prefixes <= 2 * MAX_PREFIXES_END else ["..."]),
                *self.prefixes[-MAX_PREFIXES_END:],
            ]
        else:
            print_prefixes_list = self.prefixes
        print_prefixes = ", ".join(print_prefixes_list)
        return (
            f"{self.__class__.__name__}[{self._repr_other_dims_info}]"
            f" for {self.num_prefixes} experiments [{print_prefixes}]"
        )

    def rename_prefixes(self: SelfT, prefixes: Union[str, Sequence[str]]) -> SelfT:
        """Rename prefixes as specified"""
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        if len(prefixes) != self.num_prefixes:
            raise ValueError(
                f"Given {len(prefixes)} new prefixes not matching number of"
                f" experiments in {self}"
            )
        return self._with_updated_df(self.df.assign_coords(prefix=prefixes))

    def drop_prefixes(self: SelfT, prefixes: Union[str, Collection[str]]) -> SelfT:
        if isinstance(prefixes, str):
            prefixes = {prefixes}
        else:
            prefixes = set(prefixes)
        if missing := prefixes - set(self.prefixes):
            raise ValueError(
                f"Cannot drop prefixes missing from {self}"
                f" ({missing = }, from = {self.prefixes})"
            )
        return self[[x for x in self.prefixes if x not in prefixes]]

    def __getitem__(self: SelfT, prefixes) -> SelfT:
        """Subset class to selected prefixes"""
        if isinstance(prefixes, str):
            # make sure that prefixes is a sequence
            prefixes = [prefixes]
        try:
            df = self.df.sel(prefix=prefixes)
        except KeyError as err:
            missing = set(prefixes) - set(self.prefixes)
            raise KeyError(
                "Unable to select all input prefixes"
                f" ({missing = }, from = {self.prefixes})"
            ) from err
        return self._with_updated_df(df)

    def subset_mask(self: SelfT, prefix_mask: npt.ArrayLike) -> SelfT:
        """Subset class to selected prefixes (provided as boolean mask)"""
        prefix_mask = np.array(prefix_mask, dtype=np.bool_, copy=False)
        if prefix_mask.ndim != 1:
            raise ValueError(f"prefix_mask must be 1D ({prefix_mask.shape = })")
        if len(prefix_mask) != self.num_prefixes:
            raise ValueError(
                "prefix_mask must have same length as num_prefixes"
                f" ({len(prefix_mask) = }, {self.num_prefixes = })"
            )
        return self._with_updated_df(self.df.isel(prefix=prefix_mask))

    def split_prefixes(
        self: SelfT, rng: Optional[np.random.RandomState] = None
    ) -> Tuple[SelfT, SelfT]:
        """Split class randomly into evenly sized parts"""
        if self.num_prefixes < 2:
            raise ValueError("Cannot split unless there are at least 2 experiments")
        split_mask = (rng.permutation if rng else np.random.permutation)(
            self.num_prefixes
        ) < self.num_prefixes // 2
        psi1 = self.subset_mask(split_mask)
        psi2 = self.subset_mask(~split_mask)
        return psi1, psi2

    @classmethod
    def concat(
        cls,
        *objs: SelfT,
        override_args: Optional[Sequence] = None,
        update_kwargs: Optional[Dict] = None,
    ) -> SelfT:
        """Concatenate multiple instances of class into single one

        Requires at least one object to be passed in.
        Additional args/kwargs for constructor are obtained from first object,
        but can be overriden by `override_args` and `update_kwargs`.

        Parameters
        ----------
        *objs:
            Combine the data from these input objects, using additional
            args/kwargs from the first element. At least one item must be
            specified.
        override_args: Optional[Sequence]
            Instead of using `self._non_df_args` for additional positional
            arguments, use these positional arguments
        update_kwargs: Dict
            Add or replace elements of `self._non_df_kwargs` with values in
            input dictionary for keyword arguments to constructor
        """
        if update_kwargs is None:
            update_kwargs = {}
        if not objs:
            raise ValueError("concat() requires at least one argument")
        df = xr.concat(
            [x.df for x in objs],
            "prefix",
            data_vars="minimal",
            coords="minimal",
            compat="no_conflicts",
            join="override",
            combine_attrs="drop",
        )
        return objs[0]._with_updated_df(
            df, override_args=override_args, update_kwargs=update_kwargs
        )

    def downsample(
        self: SelfT, num_prefixes: int, rng: Optional[np.random.RandomState] = None
    ) -> SelfT:
        """Get random subset with exactly `num_prefixes` prefixes

        Parameters
        ----------
        num_prefixes: int
            The number of prefixes wanted in result. Must be greater than zero.
            If `num_prefixes` greater than `self.num_prefixes`, return`self`.
        rng: np.random.RandomState, optional
            Random state used to pick which prefixes are selected.
            If not set, use numpy global random state.
        """
        if num_prefixes < 1:
            raise ValueError("Must select at least 1 prefix")
        elif num_prefixes >= self.num_prefixes:
            return self
        else:
            return self.subset_mask(
                (rng.permutation if rng else np.random.permutation)(self.num_prefixes)
                < num_prefixes
            )
