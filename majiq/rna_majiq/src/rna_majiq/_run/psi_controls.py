"""
psi_controls.py

Compute summary of input PsiCoverage with given quantiles

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_range_factory,
    resources_args,
    select_prefixes_args,
)
from rna_majiq._run._run import GenericSubcommand

DESCRIPTION = "Summarize PSI controls for repeated comparison to cases"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "controls",
        metavar="controls_summary",
        type=ResolvedPath,
        help="Path for output PsiControlsSummary file",
    )
    parser.add_argument(
        "psi",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to PsiCoverage/PsiGroup files used as controls",
    )
    select_prefixes_args(parser)
    parser.add_argument(
        "--alpha",
        metavar="A",
        # floats on [0, 1]
        type=check_range_factory(float, 0, 1, True, True),
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=nm.constants.DEFAULT_OUTLIERS_ALPHA,
        help="Threshold for control quantiles in both directions (two-sided"
        " comparison). Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha."
        " (default: %(default)s)",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.controls.exists():
        raise ValueError(
            f"Output PsiControlsSummary {args.controls} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = nm.logger.get_logger()
    log.info("Joining %d input files with PSI information", len(args.psi))
    psi = nm.PsiGroup.from_zarr_permissive(
        args.psi,
        always_convert=False,
        chunk_like_group=True,
        ec_idx_nchunks=8,
        show_progress=args.show_progress,
    )
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psi = psi[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psi = psi.drop_prefixes(args.drop_prefixes)
    log.info("Summarizing %s to %s", psi, args.controls)
    log.debug("Prefixes: %s", psi.prefixes)
    nm.PsiControlsSummary.from_psi(psi, args.alpha).to_zarr(
        args.controls, show_progress=args.show_progress
    )
    return


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args(sys_args)
    run(args)
    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


if __name__ == "__main__":
    main()
