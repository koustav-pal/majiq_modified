"""
combine.py

Combine input splicegraphs

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq._run.build_args import ir_filtering_args
from rna_majiq.logger import get_logger

DESCRIPTION = "Combine input splicegraphs into single splicegraph"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "out_sg",
        type=ResolvedPath,
        help="Path for output splicegraph",
    )
    StoreSGPaths = StoreRequiredUniqueActionFactory()
    input_args = parser.add_argument_group(
        "required splicegraph arguments (need at least one)"
    )
    input_args.add_argument(
        "--make-annotated",
        metavar="SG",
        type=ExistingResolvedPath,
        action=StoreSGPaths,
        nargs="+",
        default=list(),
        help="Input splicegraphs for which all junctions will be marked as"
        " annotated (i.e. not denovo). This helps highlight denovo junctions"
        " that were unique to non-base splicegraphs. Note that introns remain"
        " denovo in order to accommodate assignment of intronic coverage.",
    )
    input_args.add_argument(
        "--keep-denovo",
        metavar="SG",
        type=ExistingResolvedPath,
        action=StoreSGPaths,
        nargs="+",
        default=list(),
        help="Input splicegraphs for which junctions will remain marked as"
        " denovo (unless also found in inputs from --make-annotated)",
    )

    ir_filtering_args(parser)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.out_sg.exists():
        raise ValueError(
            f"Output splicegraph {args.out_sg} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    sg = nm.SpliceGraph.combine(
        args.make_annotated, args.keep_denovo, filter_introns=args.introns
    )
    log.info("Saving updated %s to %s", sg, args.out_sg)
    sg.to_zarr(args.out_sg, mode="w")
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
