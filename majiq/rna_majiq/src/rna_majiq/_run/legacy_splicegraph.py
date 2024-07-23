"""
legacy_splicegraph.py

Output splicegraph file compatible with VOILA v2

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import ExistingResolvedPath, ResolvedPath
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.logger import get_logger

DESCRIPTION = (
    "Output legacy splicegraph file (splicegraph.sql) compatible with VOILA v2"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to new-majiq splicegraph",
    )
    parser.add_argument(
        "legacy_splicegraph",
        type=ResolvedPath,
        help="Path for output legacy splicegraph file",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.legacy_splicegraph.exists():
        raise ValueError(
            f"Output splicegraph {args.legacy_splicegraph} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Loading input splicegraph from %s", args.splicegraph)
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info("Saving legacy splicegraph to %s", args.legacy_splicegraph)
    sg.to_sqlite(args.legacy_splicegraph)
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
