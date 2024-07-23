"""
heterogen.py

Analysis of groups of independent experiments

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from typing import Any, Dict, List, Optional

import numpy as np

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    check_nonnegative_factory,
    quantify_comparison_args,
    resources_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.logger import get_logger

DESCRIPTION = "Test differences in PSI for two groups of independent experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--population-quantiles",
        metavar="Q",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_POPULATION_QUANTILES,
        help="Quantiles of PSI (besides median) per group to report"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--stats",
        metavar="S",
        type=str,
        nargs="+",
        choices=set(nm.constants.STATS_AVAILABLE.keys()),
        default=nm.constants.DEFAULT_HET_USESTATS,
        help="Specify which stats to use for testing differences between groups"
        f" (available: {set(nm.constants.STATS_AVAILABLE.keys())})"
        " (default: %(default)s)",
    )
    quantify_comparison_args(parser, psi_allowed_types="PsiCoverage/PsiGroup")

    psisample_settings = parser.add_argument_group(
        "psisamples (posterior samples) testing arguments"
    )
    psisample_settings.add_argument(
        "--psisamples",
        metavar="J",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_HET_PSISAMPLES,
        help="Number of repeated samples from approximate bootstrap posterior"
        " to take, test, summarize (see --pvalue-quantiles)."
        " If set to zero, testing on approximate bootstrap posterior"
        " will only be done for posterior means (default: %(default)s).",
    )
    psisample_settings.add_argument(
        "--pvalue-quantiles",
        metavar="Q",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_PVALUE_QUANTILES,
        help="Report these quantiles of pvalues on psisamples (default: %(default)s)",
    )

    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    population_quantiles = sorted(set(np.round(args.population_quantiles, 3)))
    pvalue_quantiles = sorted(set(np.round(args.pvalue_quantiles, 3)))
    use_stats = sorted(set(args.stats))

    log = get_logger()
    if args.client:
        # count number of threads
        # workers need to have rng pool expanded to match number of available threads
        # workers need to have different random seeds
        nthreads: int = 0  # accumulate number of threads
        rng_setup = []  # futures for setting up rngs
        for worker, worker_nthreads in args.client.nthreads().items():
            if args.psisamples and args.pvalue_quantiles:
                rng_setup.append(
                    args.client.submit(
                        nm.rng_resize, worker_nthreads, workers=worker, pure=False
                    )
                )
                rng_setup.append(
                    args.client.submit(
                        nm.rng_seed,
                        args.rng_seed + nthreads,
                        workers=worker,
                        pure=False,
                    )
                )
            nthreads += worker_nthreads
        args.client.gather(rng_setup)  # wait for the rng setup to finish
    else:
        # single-machine nthreads, set up rngs
        nm.rng_resize(args.nthreads)
        nm.rng_seed(args.rng_seed)
    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__
    metadata["population_quantiles"] = population_quantiles
    metadata["stats"] = use_stats
    metadata["pvalue_quantiles"] = pvalue_quantiles
    metadata["psisamples"] = args.psisamples

    # load psi1, psi2
    log.info("Loading input groups")
    # if all PsiCoverage, do not convert yet/do not chunk like group
    # this allows downsampling/splitting to be done with more ease
    psi1 = nm.PsiGroup.from_zarr_permissive(
        args.psi1,
        always_convert=False,
        chunk_like_group=False,
        show_progress=args.show_progress,
    )
    if args.select_grp1_prefixes:
        log.debug("Selecting specified grp1 prefixes")
        psi1 = psi1[args.select_grp1_prefixes]
    elif args.drop_grp1_prefixes:
        log.debug("Dropping specified grp1 prefixes")
        psi1 = psi1.drop_prefixes(args.drop_grp1_prefixes)
    if args.psi2:
        psi2 = nm.PsiGroup.from_zarr_permissive(
            args.psi2,
            always_convert=False,
            chunk_like_group=False,
            show_progress=args.show_progress,
        )
        if args.select_grp2_prefixes:
            log.debug("Selecting specified grp2 prefixes")
            psi2 = psi2[args.select_grp2_prefixes]
        elif args.drop_grp2_prefixes:
            log.debug("Dropping specified grp2 prefixes")
            psi2 = psi2.drop_prefixes(args.drop_grp2_prefixes)
        if args.downsample2 and psi2.num_prefixes > psi1.num_prefixes:
            log.info(
                "Downsampling group 2 (%s) to have same number of experiments as"
                " group 1 (%s)",
                psi2,
                psi1,
            )
            rng = np.random.RandomState(args.rng_seed + psi2.num_prefixes)
            psi2 = psi2.downsample(psi1.num_prefixes, rng=rng)
            log.info("Group 2: %s", psi2.prefixes)
    else:
        if psi1.num_prefixes < 2:
            raise ValueError(
                f"Cannot split {psi1} unless there are at least 2 experiments"
            )
        log.info("Splitting %s randomly into two groups", psi1)
        rng = np.random.RandomState(args.rng_seed + psi1.num_prefixes)
        psi1, psi2 = psi1.split_prefixes(rng=rng)
        log.info("Group 1: %s", psi1.prefixes)
        log.info("Group 2: %s", psi2.prefixes)
    if not args.allow_prefix_overlap and (
        overlap := set(psi1.prefixes) & set(psi2.prefixes)
    ):
        raise ValueError(
            f"Group 1 ({psi1}) and 2 ({psi2}) have overlapping prefixes {overlap}"
        )
    log.debug("Group 1: %s", psi1.prefixes)
    log.debug("Group 2: %s", psi2.prefixes)

    if not psi1.events_df.equals(psi2.events_df):
        raise ValueError("Events from psi1 do not match events from psi2")
    group_sizes = {
        args.names[0]: psi1.num_prefixes,
        args.names[1]: psi2.num_prefixes,
    }
    metadata["group_sizes"] = group_sizes
    log.info(f"Comparing {args.names[0]}({psi1}) vs {args.names[1]}({psi2})")

    # dataset of quantifications I want
    log.info("Preparing distribution parameters for testing")
    heterogen = nm.Heterogen(
        psi1,
        psi2,
        min_experiments_f=args.min_experiments,
        name1=args.names[0],
        name2=args.names[1],
        show_progress=args.show_progress,
    )
    log.info(f"Performing quantification on means and {args.psisamples} psisamples")
    heterogen_voila = heterogen.dataset(
        pvalue_quantiles=pvalue_quantiles,
        use_stats=use_stats,
        psisamples=args.psisamples,
    )

    if args.output_voila:
        log.info("Saving quantifications for VOILA to %s", args.output_voila)
        heterogen_voila.to_zarr(
            args.output_voila,
            show_progress=args.show_progress,
        )
        heterogen_voila = nm.HeterogenDataset.from_zarr(args.output_voila)

    sg: Optional[nm.SpliceGraph] = None
    annotated: Optional[nm.SpliceGraph] = None
    if args.splicegraph:
        log.debug("Loading splicegraph from %s", args.splicegraph)
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)
        if args.annotated:
            log.debug("Loading splicegraph (annotated) from %s", args.annotated)
            annotated = nm.SpliceGraph.from_zarr(args.annotated, genes=sg.genes)

    log.info("Summarizing population quantiles per group")
    df = heterogen_voila.to_dataframe(
        sg=sg,
        annotated=annotated,
        population_quantiles=population_quantiles,
        show_progress=args.show_progress,
    )

    try:
        output_name = args.output_tsv.name
    except AttributeError:
        output_name = args.output_tsv
    log.info(f"Writing metadata to {output_name}")
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info(f"Writing table to {output_name}")
    (
        df
        # any column with pvalue in it needs to be manually formatted
        .pipe(
            lambda df: df.assign(
                **{
                    col: df[col].apply(lambda x: f"{x:.3e}")
                    for col in df.columns
                    if "pvalue" in col
                }
            )
        )
        # other numeric columns need at most 4 digits precision (psi)
        .round(4)
        # save result in TSV format
        .to_csv(args.output_tsv, sep="\t", index=not args.splicegraph)
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
