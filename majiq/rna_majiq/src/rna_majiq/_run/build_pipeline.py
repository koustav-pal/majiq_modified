"""
build_pipeline.py

Run commands of majiq-build in sequence:
+ Load annotated splicegraph (majiq-build gff3)
+ Load coverage from BAM files (majiq-build sj)
+ Update annotated splicegraph (majiq-build update)
"""

import argparse
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Any, Callable, List, Optional

import pandas as pd

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    ResolvedPath,
    resources_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq._run.build_args import (
    SJGroupsT,
    build_threshold_args,
    build_type_args,
    do_build_and_simplify,
    enable_simplify_args,
    gff3_parsing_args,
    gff3_parsing_pipeline,
    ir_filtering_args,
    simplifier_threshold_args,
    sj_strandness_args,
)
from rna_majiq.experiments import bam_experiment_name
from rna_majiq.logger import get_logger

DESCRIPTION = (
    "majiq-build pipeline to build splicegraph"
    " from annotations and RNA-seq experiments"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "gff3",
        type=ExistingResolvedPath,
        help="Path to GFF3 file (uncompressed or gzipped) with annotated"
        " transcripts for initial splicegraph",
    )
    parser.add_argument(
        "experiments_tsv",
        type=ExistingResolvedPath,
        help="Path to TSV file indicating experiments from one or more build"
        " groups."
        " Required column 'path', indicating path to input file for each unique"
        " experiment."
        " Column `group` indicates group experiment belongs to; `path`"
        " indicates the path to the input file used for the experiment."
        " Optional columns 'group' and 'strandness'."
        " Column `group` indicates group experiment belongs to, by default all"
        " experiments are in a single build group."
        " Column `strandness` overrides value from --strandness."
        " It is assumed that previously loaded SJExperiment have extension"
        " `sj`, `zarr`, or `zip` (and alignments do not).",
    )
    parser.add_argument(
        "output_dir",
        type=ResolvedPath,
        help="Path for output directory with files produced by the pipeline",
    )
    build_type_args(parser, allow_simplify_only=False)
    ir_filtering_args(parser)
    enable_simplify_args(parser)
    build_threshold_args(parser)
    simplifier_threshold_args(parser, prefix="simplify-")
    sj_strandness_args(parser)
    gff3_parsing_args(parser, prefix="gff3-")
    resources_args(parser, use_dask=False)
    return


def load_experiments_tsv(
    path: Path,
    global_strandness: str,
    output_dir: Path,
    path_resolver: Callable[[Any], Path] = NewResolvedPath,
) -> pd.DataFrame:
    """Load table of experiments, validating paths, setting defaults"""
    log = get_logger()
    df_experiments = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "path": str,
            "group": str,
            "strandness": str,
            "is_sj": bool,
        },
    )

    # validate experiments_tsv
    if "path" not in df_experiments.columns:
        raise ValueError("Experiments table does not have 'path' column")
    elif len(set(df_experiments["path"])) < len(df_experiments):
        nonunique_paths = (
            df_experiments["path"]
            .value_counts()
            .pipe(lambda x: x[x > 1])
            .index.tolist()
        )
        raise ValueError(f"Input paths are not unique ({nonunique_paths = })")
    else:
        # verifies that all input paths exist
        df_experiments["path"] = df_experiments["path"].apply(ExistingResolvedPath)

    if "group" not in df_experiments.columns:
        log.debug("No `group` column ; all experiments will be in one build group")
        df_experiments["group"] = "experiments"  # all are same group, trivial name
    elif df_experiments["group"].isna().any():
        raise ValueError("Specified group column may not have missing values")

    if "strandness" not in df_experiments.columns:
        df_experiments["strandness"] = global_strandness
    else:
        df_experiments.fillna({"strandness": global_strandness}, inplace=True)
        # all should be uppercase
        df_experiments["strandness"] = df_experiments["strandness"].str.upper()
        valid_strandness = {"AUTO", *nm.ExperimentStrandness.__entries.keys()}
        if invalid_strandness := set(df_experiments["strandness"]) - valid_strandness:
            raise ValueError(
                "Invalid values for strandness in experiments_tsv"
                f" ({invalid_strandness = }, {valid_strandness = })"
            )

    df_experiments["is_sj"] = df_experiments["path"].apply(
        lambda x: x.suffix.upper() in {".SJ", ".ZIP", ".ZARR"}
    )
    df_experiments["sj_path"] = [
        # verifies the SJ files created from BAM do not already exist
        path if is_sj else path_resolver(output_dir / f"{bam_experiment_name(path)}.sj")
        for path, is_sj in zip(df_experiments["path"], df_experiments["is_sj"])
    ]

    return df_experiments


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info("Loading experiment information from %s", args.experiments_tsv)

    # path resolver
    path_resolver = ResolvedPath if args.overwrite else NewResolvedPath

    # load args.experiments_tsv
    df_experiments = load_experiments_tsv(
        args.experiments_tsv,
        args.strandness,
        args.output_dir,
        path_resolver=path_resolver,
    )
    # path for output splicegraph
    output_splicegraph = path_resolver(args.output_dir / "splicegraph.zarr")

    log.info("Preparing directory for output files at %s", args.output_dir)
    args.output_dir.mkdir(exist_ok=True)

    # load base splicegraph
    log.info("(`majiq-build gff3`)")
    sg: nm.SpliceGraph = gff3_parsing_pipeline(
        args.gff3,
        features_default=args.gff3_features_default,
        features_tsv=args.gff3_features_tsv,
        types_genes=args.gff3_types_genes,
        types_transcripts=args.gff3_types_transcripts,
        types_exons=args.gff3_types_exons,
        types_silent=args.gff3_types_silent,
        types_hard_skip=args.gff3_types_hard_skip,
        save_annotated=not args.skip_saving_annotated_transcripts,
    )
    log.info("Loaded annotated %s", sg)

    # get SJ files that we don't already have
    log.info("(`majiq-build sj`)")
    sj: Optional[nm.SJExperiment] = None
    df_experiments_bam = df_experiments.loc[~df_experiments["is_sj"]]
    for idx, (_, experiment_info) in enumerate(df_experiments_bam.iterrows(), 1):
        sj = nm.SJExperiment.from_bam(
            experiment_info["path"],
            sg,
            experiment_info["strandness"],
            nthreads=args.nthreads,
            auto_minreads=args.auto_minreads,
            auto_minjunctions=args.auto_minjunctions,
            auto_mediantolerance=args.auto_mediantolerance,
        )
        log.info(
            "Saving %s to %s (%d / %d)",
            sj,
            experiment_info["sj_path"],
            idx,
            len(df_experiments_bam),
        )
        sj.to_zarr(experiment_info["sj_path"], mode="w")
    sj = None  # unload SJ coverage from memory

    # create pool for multithreading
    p = ThreadPool(args.nthreads)

    log.info("(`majiq-build update`)")
    experiments: SJGroupsT = {
        group: sj_paths.tolist()
        for group, sj_paths in df_experiments.groupby("group")["sj_path"]
    }
    sg = do_build_and_simplify(
        sg,
        experiments,
        build=args.build,
        minreads=args.minreads,
        mindenovo=args.mindenovo,
        minpos=args.minpos,
        max_pctbins=args.max_pctbins,
        junction_acceptance_probability=args.junction_acceptance_probability,
        intron_acceptance_probability=args.intron_acceptance_probability,
        min_experiments=args.min_experiments,
        introns_type=args.introns,
        simplify=bool(args.simplify),
        reset_simplify=True,
        simplify_minpsi=args.simplify_minpsi,
        simplify_minreads_annotated=args.simplify_minreads_annotated,
        simplify_minreads_denovo=args.simplify_minreads_denovo,
        simplify_minreads_ir=args.simplify_minreads_ir,
        simplify_min_experiments=args.simplify_min_experiments,
        imap_unordered_fn=p.imap_unordered,
    )
    log.info("Saving updated %s to %s", sg, output_splicegraph)
    sg.to_zarr(output_splicegraph, mode="w")

    p.close()
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
