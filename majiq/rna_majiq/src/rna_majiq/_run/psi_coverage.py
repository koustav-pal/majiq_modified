"""
psi_coverage.py

Get coverage over LSVs for a specific experiment

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
    prefixes_args,
    resources_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq._run.build_args import lsv_coverage_args, quantifiability_threshold_args
from rna_majiq.logger import get_logger

DESCRIPTION = "Prepare raw and bootstrapped coverage at LSVs for quantification"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments for parser"""
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph to define LSVs",
    )
    parser.add_argument(
        "psi_coverage",
        type=ResolvedPath,
        help="Path for output psi-coverage file",
    )
    parser.add_argument(
        "sj",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to SJ coverage files for experiments",
    )
    prefixes_args(parser)
    quantifiability_threshold_args(parser)
    lsv_coverage_args(parser)
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    if not args.overwrite and args.psi_coverage.exists():
        raise ValueError(
            f"Output PsiCoverage {args.psi_coverage} already exists"
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
    log.debug(f"Defining events for coverage ({args.select_lsvs})")
    lsvs = sg.exon_connections.lsvs(args.select_lsvs)
    if args.ignore_from is not None:
        log.info(f"Ignoring events also found in {args.ignore_from}")
        lsvs = lsvs[
            lsvs.unique_events_mask(
                nm.SpliceGraph.from_zarr(
                    args.ignore_from, genes=sg.genes
                ).exon_connections.lsvs(args.select_lsvs)
            ).unique_events_mask
        ]
    log.info("Assessing coverage for %s on %s", lsvs, sg)

    # if ZIP output, create unzipped intermediate for batched creation of zarr
    if str(args.psi_coverage).endswith(".zip"):
        batch_output = args.work_directory / f"{args.psi_coverage.name}.unzipped"
        log.info("Writing unzipped PsiCoverage to working path %s", batch_output)
    else:
        batch_output = args.psi_coverage

    nm.rng_resize(args.nthreads)
    p: Optional[ThreadPool] = None
    imap_unordered_fn: Callable = map
    if len(args.sj) != 1:
        p = ThreadPool(args.nthreads)
        imap_unordered_fn = p.imap_unordered

    nm.PsiCoverage.convert_sj_batch(
        args.sj,
        lsvs,
        batch_output,
        minreads=args.quantify_minreads,
        minbins=args.quantify_minbins,
        num_bootstraps=args.num_bootstraps,
        pvalue_threshold=args.stack_pvalue_threshold,
        imap_unordered_fn=imap_unordered_fn,
        prefixes=args.prefixes,
    )

    if p:
        p.close()

    if str(args.psi_coverage).endswith(".zip"):
        # make zip file output
        log.info("Zipping PsiCoverage to final path %s", args.psi_coverage)
        make_archive(str(args.psi_coverage)[:-4], "zip", batch_output)

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
