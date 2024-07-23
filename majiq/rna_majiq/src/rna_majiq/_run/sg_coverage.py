"""
sg_coverage.py

Get coverage over splicegraph for a specific experiment

Author: Joseph K Aicher
"""

import argparse
from multiprocessing.pool import ThreadPool
from shutil import make_archive
from typing import Callable, List, Optional

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
    chunks_args,
    prefixes_args,
    resources_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.logger import get_logger

DESCRIPTION = (
    "Get raw coverage for input experiment at splicegraph introns and junctions"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph with introns/junctions",
    )
    parser.add_argument(
        "sg_coverage",
        type=ResolvedPath,
        help="Path for output coverage over introns/junctions",
    )
    parser.add_argument(
        "sj",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to SJ coverage for experiments",
    )
    prefixes_args(parser)
    chunks_args(parser, nm.constants.NC_SGREADS_CHUNKS)
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.sg_coverage.exists():
        raise ValueError(
            f"Output SpliceGraphReads {args.sg_coverage} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    if args.prefixes is not None and len(args.prefixes) != len(args.sj):
        raise ValueError(
            "Using --prefixes requires equal number of prefixes and SJ files"
            f" ({len(args.prefixes) = }, {len(args.sj) = })"
        )
    log = get_logger()
    log.info(f"Loading input splicegraph from {args.splicegraph}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)

    # if ZIP output, create unzipped intermediate for batched creation of zarr
    if str(args.sg_coverage).endswith(".zip"):
        batch_output = args.work_directory / f"{args.sg_coverage.name}.unzipped"
        log.info("Writing unzipped SpliceGraphReads to working path %s", batch_output)
    else:
        batch_output = args.sg_coverage

    p: Optional[ThreadPool] = None
    imap_unordered_fn: Callable = map
    if len(args.sj) != 1:
        p = ThreadPool(args.nthreads)
        imap_unordered_fn = p.imap_unordered

    nm.SpliceGraphReads.convert_sj_batch(
        args.sj,
        sg.introns,
        sg.junctions,
        batch_output,
        chunksize=args.chunksize,
        attrs=dict(sg=str(args.splicegraph)),
        imap_unordered_fn=imap_unordered_fn,
        prefixes=args.prefixes,
    )

    if p:
        p.close()

    if str(args.sg_coverage).endswith(".zip"):
        # make zip file output
        log.info("Zipping SpliceGraphReads to final path %s", args.sg_coverage)
        make_archive(str(args.sg_coverage)[:-4], "zip", batch_output)

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
