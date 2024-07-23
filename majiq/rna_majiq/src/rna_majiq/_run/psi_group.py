"""
psi_group.py

Group PsiCoverage over prefixes for faster computation with HET, etc.
Also, summarize information from bootstrap replicates that is only necessary
for MAJIQ deltapsi.

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
    resources_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.logger import get_logger

DESCRIPTION = (
    "Summarize posteriors for efficient computation over prefixes (HET, Controls)"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "output",
        type=ResolvedPath,
        help="Path for output PsiGroup file",
    )
    parser.add_argument(
        "psi",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to PsiCoverage/PsiGroup files that should be made part of the group",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.output.exists():
        raise ValueError(
            f"Output PsiGroup {args.output} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Joining %d input files with PSI information", len(args.psi))
    psi = nm.PsiGroup.from_zarr_permissive(
        args.psi,
        always_convert=False,
        # if all PsiCoverage, chunk as PsiCoverage for efficient conversion
        chunk_like_group=False,
        show_progress=args.show_progress,
    )
    log.debug("Prefixes: %s", psi.prefixes)
    if isinstance(psi, nm.PsiCoverage):
        log.info("Converting %s to PsiGroup; saving to %s", psi, args.output)
        psi = psi.group(save_zarr=args.output, show_progress=args.show_progress)
    else:
        log.info("Saving %s to %s", psi, args.output)
        psi.to_zarr(args.output)
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
