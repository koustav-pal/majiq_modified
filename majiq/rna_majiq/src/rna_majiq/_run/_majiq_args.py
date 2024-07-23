"""
majiq_args.py

Helper classes/functions for checking valid arguments directly at command-line
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Union

import rna_majiq as nm


def ResolvedPath(x: Any) -> Path:
    try:
        p = Path(x).resolve()
    except TypeError as err:
        raise argparse.ArgumentTypeError(f"Input {x} is not pathlike") from err
    return p


def ExistingResolvedPath(x: Any) -> Path:
    p = ResolvedPath(x)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"{p} does not exist")
    return p


def NewResolvedPath(x: Any) -> Path:
    p = ResolvedPath(x)
    if p.exists():
        raise argparse.ArgumentTypeError(f"{p} already exists")
    return p


def ExistingResolvedDirectory(x: Any) -> Path:
    p = ResolvedPath(x)
    if not p.is_dir():
        raise argparse.ArgumentTypeError(f"{p} is not a directory")
    return p


def check_characters_factory(
    valid_characters: Set[str],
    valid_description: Optional[str] = None,
) -> Callable[[Any], str]:
    """Create function to raise argparse type error or return string

    Create function to raise argparse.ArgumentTypeError if input contains
    invalid characters, otherwise, return str(x)

    Parameters
    ----------
    valid_characters: Set[str]
        characters that are acceptable for this restricted string
    valid_description: Optional[str]
        description of acceptable characters in input
    """

    def check_characters(x: Any) -> str:
        x = str(x)
        set_x = set(x)
        if not set_x.issubset(valid_characters):
            invalid_chars = set_x.difference(valid_characters)
            prefix = (
                f"may only contain {valid_description}: " if valid_description else ""
            )
            raise argparse.ArgumentTypeError(
                f"{prefix}{x} contains invalid characters {invalid_chars}"
            )
        return x

    return check_characters


def check_range_factory(
    cast_fn: Callable[[Any], Union[int, float]],
    minval: Union[int, float] = float("-inf"),
    maxval: Union[int, float] = float("inf"),
    mininclude: bool = False,
    maxinclude: bool = False,
):
    """Return argparse type function to casted type in specified range"""
    not_too_small = (lambda x: x >= minval) if mininclude else (lambda x: x > minval)
    not_too_big = (lambda x: x <= maxval) if maxinclude else (lambda x: x < maxval)
    invalid_message = (
        "parameter must be in range"
        f" {'[' if mininclude else '('}{minval}, {maxval}{']' if maxinclude else ')'}"
    )

    def return_function(value):
        ivalue = cast_fn(value)
        if not_too_small(ivalue) and not_too_big(ivalue):
            return ivalue
        else:
            raise argparse.ArgumentTypeError(
                f"{value} outside of valid range; {invalid_message}"
            )

    return return_function


def check_nonnegative_factory(
    cast_fn: Callable[[Any], Union[int, float]], reject_zero: bool
):
    """Returns argparse type function to casted type, rejecting negative values"""
    return check_range_factory(cast_fn, minval=0, mininclude=not reject_zero)


def StoreRequiredUniqueActionFactory(append: bool = False):
    """Return class that tracks and prevents repeated values

    Parameters
    ----------
    append: bool
        If append is True, treat as an 'append' action, so that input values
        are appended to destination list
    """

    class StoreRequiredUniqueAction(argparse.Action):
        # class variable with shared values known to all instances
        shared_values: Dict[Any, Set] = {}

        @staticmethod
        def repeated_values(values) -> Set:
            """Returns repeated values if input is list (empty set if not list)"""
            if isinstance(values, list):
                unique_values = set(values)
                if len(values) != len(unique_values):
                    repeated_list = values.copy()
                    for x in unique_values:
                        repeated_list.remove(x)
                    return set(repeated_list)
            return set()

        @classmethod
        def overlapping_values(
            cls, dest, values, update: bool = True
        ) -> Dict[Any, Set]:
            """Returns values overlapping with previously processed arguments"""
            unique_values = set(values)
            overlapping_values = {
                other_key: other_values & unique_values
                for other_key, other_values in cls.shared_values.items()
                if other_values & unique_values
            }
            if update:
                try:
                    cls.shared_values[dest] |= unique_values
                except KeyError:
                    cls.shared_values[dest] = unique_values
            return overlapping_values

        def _check_unique_and_update(self, values: List) -> None:
            """Checks if values are unique/overlapping and update overlapping"""
            # check for repeated values
            if repeated_values := self.repeated_values(values):
                raise argparse.ArgumentError(
                    self, f"Values to {self.dest} not unique ({repeated_values = })"
                )
            # check if overlapping previously processed values, update with new values
            if overlapping_values := self.overlapping_values(self.dest, values):
                raise argparse.ArgumentError(
                    self,
                    f"Values to {self.dest} overlap previously processed arguments"
                    f" ({overlapping_values = })",
                )
            return

        def __call__(self, parser, namespace, values, option_string=None):
            """Check that values are unique if list, not in shared_values"""
            if isinstance(values, list):
                self._check_unique_and_update(values)
            else:
                self._check_unique_and_update([values])
            # update namespace
            if append:
                try:
                    getattr(namespace, self.dest).append(values)
                except AttributeError:
                    setattr(namespace, self.dest, [values])
            else:
                setattr(namespace, self.dest, values)
            return

    return StoreRequiredUniqueAction


def _quantify_shared_args(parser: argparse.ArgumentParser) -> None:
    """shared arguments for quantify_?(no)comparison_args"""
    parser.add_argument(
        "--splicegraph",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="If specified, annotate quantifications with splicegraph information",
    )
    parser.add_argument(
        "--annotated",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="If specified, identify novel events/exons/introns relative to"
        " this splicegraph (vs --splicegraph)."
        " Ignored if --splicegraph not used.",
    )
    parser.add_argument(
        "--output-tsv",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )


def quantify_nocomparison_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to PsiCoverage files to quantify",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        default=None,
        help="If specified, treat samples as replicates and quantify combined"
        " coverage for events that passed in at least min experiments"
        " (proportion of experiments if < 1)."
        " Default is to quantify each sample independently.",
    )
    select_prefixes_args(parser)
    _quantify_shared_args(parser)
    return


def quantify_comparison_args(
    parser: argparse.ArgumentParser,
    psi_allowed_types: str = "PsiCoverage",
) -> None:
    comparison_req = parser.add_argument_group("quantification group arguments")
    comparison_req.add_argument(
        "-psi1",
        "-grp1",
        "-g1",
        dest="psi1",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        required=True,
        help=f"Paths to {psi_allowed_types} files for experiments contributing"
        " to psi1",
    )
    comparison_req.add_argument(
        "-psi2",
        "-grp2",
        "-g2",
        dest="psi2",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        required=False,
        default=None,
        help=f"Paths to {psi_allowed_types} files for experiments contributing"
        " to psi2."
        " If not specified, the experiments contributing to psi1 will be"
        " split into two groups for the comparison.",
    )
    StoreGroupNames = StoreRequiredUniqueActionFactory()
    check_group_chars = check_characters_factory(
        nm.constants.ALLOWED_GROUP_NAME_CHARS, "alphanumeric or underscore characters"
    )
    comparison_req.add_argument(
        "-n",
        "--names",
        dest="names",
        nargs=2,
        metavar=("NAME1", "NAME2"),
        default=["grp1", "grp2"],
        action=StoreGroupNames,
        type=check_group_chars,
        help="The names that identify the groups being compared"
        " (default: %(default)s).",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. This specifies the fraction"
        " (value < 1) or absolute number (value >= 1) of experiments that must"
        " pass individually in each group for an LSV to be quantified"
        " (default: %(default)s)",
    )
    _quantify_shared_args(parser)
    parser.add_argument(
        "--output-voila",
        metavar="ZARR",
        type=NewResolvedPath,
        default=None,
        help="Path for output voila view file (default: none)",
    )
    select_prefixes_args(parser, "grp1")
    select_prefixes_args(parser, "grp2")
    parser.add_argument(
        "--downsample2",
        action="store_true",
        default=False,
        help="If psi2 has more experiments than psi2, a random subset of"
        " experiments from psi2 will be selected to perform the comparison"
        " on groups of equal size (default: no downsampling).",
    )
    parser.add_argument(
        "--rng-seed",
        metavar="S",
        type=int,
        default=22022022,
        help="Set seed for random number generators (default: %(default)s)",
    )
    parser.add_argument(
        "--allow-prefix-overlap",
        default=False,
        action="store_true",
        help="Allow groups compared to have overlapping prefixes"
        " (default: %(default)s)",
    )
    return


def chunks_args(parser: argparse.ArgumentParser, chunksize: Optional[int]) -> None:
    chunks = parser.add_argument_group("output chunks arguments")
    if chunksize is None:
        chunks.add_argument(
            "--chunksize",
            metavar="CHUNKS",
            type=check_nonnegative_factory(int, True),
            default=chunksize,
            help="If specified, chunk output this many connections at a time."
            " Use to enforce matching chunksizes over multiple analyses."
            " By default, let MAJIQ pick chunksize automatically.",
        )
    else:
        chunks.add_argument(
            "--chunksize",
            metavar="CHUNKS",
            type=check_nonnegative_factory(int, True),
            default=chunksize,
            help="Chunk output this many introns/junctions at a time"
            " (default: %(default)s)",
        )


def select_prefixes_args(parser: argparse.ArgumentParser, group_name: str = "") -> None:
    """add exclusive arguments to select/drop prefixes

    Add exlusive arguments to select/drop prefixes. If not group_name,
    flags `--select-prefixes / --drop-prefixes`.
    If group_name specified, `--{select,drop}-[group]-prefixes`
    """
    group_description_suffix = f" ({group_name})" if group_name else ""
    select_group = parser.add_argument_group(
        f"select prefixes arguments{group_description_suffix}"
    )
    group_arg = f"-{group_name}-" if group_name else "-"
    select_group_ex = select_group.add_mutually_exclusive_group()
    select_group_ex.add_argument(
        f"--select{group_arg}prefixes",
        metavar="PREFIX",
        type=str,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="If provided, select subset of prefixes for"
        f" analysis{group_description_suffix} (default: use all prefixes)",
    )
    select_group_ex.add_argument(
        f"--drop{group_arg}prefixes",
        metavar="PREFIX",
        type=str,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="If provided, drop subset of prefixes for"
        f" analysis{group_description_suffix} (default: use all prefixes)",
    )
    return


def prefixes_args(parser: argparse.ArgumentParser) -> None:
    """Set prefixes matching input SJ files (instead of what's stored in SJ file)"""
    parser.add_argument(
        "--prefixes",
        type=str,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="If provided, names that should be used as `prefixes`"
        " (or experiment names) for each input SJ file."
        " By default, use prefix from original BAM file name stored in each"
        " SJ file.",
    )
    return


def resources_args(parser: argparse.ArgumentParser, use_dask: bool = False) -> None:
    resources = parser.add_argument_group("threads/resources arguments")
    if use_dask:
        resources.add_argument(
            "--nthreads",
            "-j",
            metavar="N",
            type=check_nonnegative_factory(int, True),
            default=nm.constants.DEFAULT_QUANTIFY_NTHREADS,
            help="Number of threads to perform work in chunks for Dask scheduler"
            " (default: %(default)s)",
        )
    else:
        resources.add_argument(
            "--nthreads",
            "-j",
            metavar="N",
            type=check_nonnegative_factory(int, True),
            default=nm.constants.DEFAULT_BAM_NTHREADS,
            help="Number of threads used for simultaneous processing of multiple"
            " input files (default: %(default)s)",
        )
    # next argument will be working directory
    work_dir_ex = resources.add_mutually_exclusive_group()
    work_dir_ex.add_argument(
        "--work-directory",
        metavar="path",
        type=ExistingResolvedDirectory,
        default=None,
        help="If specified, retain any intermediate files in specified directory",
    )
    work_dir_ex.add_argument(
        "--tmp-work-directory-parent",
        metavar="path",
        type=ExistingResolvedDirectory,
        default=None,
        help="Parent of temporary directory for intermediate files (default: chosen by Python `tempfile`)",
    )
    if use_dask:
        progress_bar = resources.add_mutually_exclusive_group()
        progress_bar.add_argument(
            "--show-progress",
            dest="show_progress",
            action="store_true",
            default=True,
            help="Show progress bar for Dask computations"
            " (default: show_progress=%(default)s)",
        )
        progress_bar.add_argument(
            "--disable-progress",
            dest="show_progress",
            action="store_false",
            default=True,
            help="Disable progress bar for Dask computations"
            " (default: show_progress=%(default)s)",
        )
        scheduler_ex = resources.add_mutually_exclusive_group(required=False)
        scheduler_ex.add_argument(
            "--scheduler-address",
            metavar="url",
            type=str,
            default=None,
            help="Use resources from existing Dask cluster at specified URL."
            " This ignores nthreads."
            " Default is to create a new local cluster with those resources"
            " for the duration of the script."
            " If used, expects cluster to remain running and accessible while"
            " the script is running.",
        )
        scheduler_ex.add_argument(
            "--scheduler-file",
            metavar="json",
            type=ExistingResolvedPath,
            default=None,
            help="Use resources from existing Dask cluster as specified in"
            " Dask scheduler file (JSON-formatted)."
            " This ignores nthreads."
            " Default is to create a new local cluster with those resources"
            " for the duration of the script."
            " If used, expects cluster to remain running and accessible while"
            " the script is running.",
        )
        parser.set_defaults(use_dask=True)
