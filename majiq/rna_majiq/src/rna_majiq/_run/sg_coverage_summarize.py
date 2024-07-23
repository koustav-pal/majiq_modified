"""
sg_coverage_summmarize.py

Summarize over sg_coverage

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
    chunks_args,
    resources_args,
    select_prefixes_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.experiments import bam_experiment_name
from rna_majiq.logger import get_logger

DESCRIPTION = "Summarize splicegraph coverage from multiple experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "summary",
        type=ResolvedPath,
        help="Path for resulting summary splicegraph reads"
        " (summary experiment name will be inferred from the path)",
    )
    parser.add_argument(
        "sg_coverage",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path for input coverage over introns/junctions",
    )
    select_prefixes_args(parser)
    parser.add_argument(
        "--reduction",
        type=str,
        choices=sorted({"sum", "mean", "median", "max"}),
        default="median",
        help="Summary function applied to coverage from input experiments (default: %(default)s)",
    )
    parser.add_argument(
        "--new-prefix",
        type=str,
        default=None,
        help="Set prefix for summarized coverage."
        " (default: use prefix of output file name)",
    )
    chunks_args(parser, nm.constants.NC_SGREADS_CHUNKS)
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.summary.exists():
        raise ValueError(
            f"Output SpliceGraphReads {args.summary} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Joining %d input splicegraph coverage files", len(args.sg_coverage))
    coverage = nm.SpliceGraphReads.from_zarr(sorted(set(args.sg_coverage)))
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        coverage = coverage[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        coverage = coverage.drop_prefixes(args.drop_prefixes)
    log.info(
        "Summarizing %s by taking the %s over all prefixes", coverage, args.reduction
    )
    log.debug("Prefixes: %s", coverage.prefixes)
    coverage = coverage.summarize(
        args.new_prefix or bam_experiment_name(args.summary), reduction=args.reduction
    )
    log.info("Saving %s to %s", coverage, args.summary)
    coverage.to_zarr(
        args.summary,
        mode="w",
        chunksize=args.chunksize,
        show_progress=args.show_progress,
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
