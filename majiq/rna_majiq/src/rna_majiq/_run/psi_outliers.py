"""
psi_outliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from typing import Any, Dict, List, Optional, Union

import pandas as pd

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    check_range_factory,
    resources_args,
    select_prefixes_args,
)
from rna_majiq._run._run import GenericSubcommand

DESCRIPTION = "Identify outliers from PSI cases (N ~ 1) vs PSI controls (N >> 1)"


def add_args(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("required analysis pass arguments")
    required.add_argument(
        "--pass",
        dest="pass_args_list",
        metavar=("splicegraph", "controls", "cases"),
        nargs=3,
        type=ExistingResolvedPath,
        action="append",
        required=True,
        help="Specify an analysis pass of SpliceGraph, PsiControlsSummary,"
        " and PsiCoverage|PsiGroup for cases for outliers analysis."
        " Use this flag multiple times to merge analysis over multiple"
        " build/quantification passes."
        " Final results will be relative to the last specified pass."
        " Each pass' PsiControlsSummary and PsiCoverage|PsiGroup must share events"
        " defined on the pass' SpliceGraph."
        " Each SpliceGraph must share genes."
        " Controls are required to share the same experiments."
        " The cases from each pass must include the same experiments as the"
        " cases from the last specified pass.",
    )
    parser.add_argument(
        "--abs-dpsi-threshold",
        type=check_range_factory(float, minval=0, mininclude=True),
        default=nm.constants.DEFAULT_OUTLIERS_DPSI_THRESHOLD,
        help="Minimum absolute difference between median PSI in controls and"
        " posterior mean PSI in case for an event connection to be considered"
        " for outlier status. Set to 0.0 to consider all connections regardless"
        " of absolute difference in PSI. (default: %(default)s)",
    )
    annotated = parser.add_argument_group("annotated features arguments")
    annotated_ex = annotated.add_mutually_exclusive_group(required=False)
    annotated_ex.add_argument(
        "--annotated",
        dest="annotated",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="Identify novel events/exons/introns/junctions relative to this"
        " splicegraph."
        " Default: use the splicegraph from the second-last analysis pass"
        " if more than one pass was specified (otherwise '--no-annotated').",
    )
    annotated_ex.add_argument(
        "--no-annotated",
        dest="annotated",
        action="store_false",
        default=None,
        help="Identify novel exons/introns/junctions using definitions from"
        " last specified SpliceGraph."
        " Do not identify novel features relative to another splicegraph."
        " Default: use the splicegraph from the second-last analysis pass"
        " if more than one pass was specified",
    )
    tsv = parser.add_argument_group("output TSV locations arguments")
    tsv.add_argument(
        "--output-tsv",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )
    tsv.add_argument(
        "--events-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing event (default: None)",
    )
    tsv.add_argument(
        "--genes-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing gene (default: None)",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_range_factory(float, minval=0, mininclude=False),
        default=nm.constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
        help="The fraction (value < 1) or absolute number (value >= 1) of"
        " experiments that must pass individually for an event in the controls"
        " for a comparison to be made (default: %(default)s)",
    )
    parser.add_argument(
        "--select-alpha",
        "--select-alpha-controls",
        metavar="A",
        type=check_range_factory(float, 0, 1, True, True),
        default=None,
        help="Select single threshold on controls quantiles in both"
        " directions (two-sided comparison)."
        " Quantiles are 0.5*alpha and 1.0-0.5*alpha."
        " Selected threshold must have been used when creating controls"
        " (default: use all values of alpha available in controls)",
    )
    parser.add_argument(
        "--select-alpha-case",
        metavar="AC",
        type=check_range_factory(float, 0, 1, True, True),
        default=None,
        help="Select single threshold on case quantiles in both"
        " directions (two-sided comparison)."
        " Quantiles are 0.5*alpha and 1.0-0.5*alpha."
        " (default: use the same values as for controls)",
    )
    select_prefixes_args(parser, "cases")
    parser.add_argument(
        "--allow-prefix-overlap",
        default=False,
        action="store_true",
        help="Allow groups compared to have overlapping prefixes"
        " (default: %(default)s)",
    )
    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    log = nm.logger.get_logger()

    # args.pass_args_list: List[Tuple[sg, controls, cases]]
    sg_paths, controls_paths, cases_paths = zip(*args.pass_args_list)

    log.debug("Loading input splicegraphs")
    use_genes = nm.Genes.from_zarr(sg_paths[-1])
    sg_list: List[nm.SpliceGraph] = [
        nm.SpliceGraph.from_zarr(sg_path, genes=use_genes) for sg_path in sg_paths
    ]
    log.info(
        "Outliers analysis on %d passes with splicegraphs %s", len(sg_list), sg_list
    )
    # main splicegraph is always splicegraph from last pass
    sg_main: nm.SpliceGraph = sg_list[-1]
    annotated: Optional[nm.SpliceGraph]
    if args.annotated:
        try:
            annotated = sg_list[sg_paths.index(args.annotated)]
        except ValueError:  # args.annotated not from one of the input passes
            annotated = nm.SpliceGraph.from_zarr(args.annotated, genes=use_genes)
        log.info(
            "Novel features in %s (%s) will be defined relative to %s (%s)",
            sg_paths[-1],
            sg_main,
            args.annotated,
            annotated,
        )
    elif args.annotated is None and len(sg_list) > 1:
        annotated = sg_list[-2]
        log.info(
            "Novel features in %s (%s) will be defined relative to %s (%s)",
            sg_paths[-1],
            sg_main,
            sg_paths[-2],
            annotated,
        )
    else:
        # only one pass OR explicit --no-annotated -> args.annotated == False
        annotated = None
        log.info(
            "Novel features in %s (%s) will be defined as noted in the splicegraph",
            sg_paths[-1],
            sg_main,
        )
    log.debug("Loading input controls")
    controls_list: List[nm.PsiControlsSummary] = [
        nm.PsiControlsSummary.from_zarr(controls_path)
        for controls_path in controls_paths
    ]
    log.debug("Loading input cases")
    cases_list: List[Union[nm.PsiCoverage, nm.PsiGroup]] = [
        nm.PsiGroup.from_zarr_permissive(
            cases_path,
            always_convert=False,
            chunk_like_group=False,
            show_progress=False,
        )
        for cases_path in cases_paths
    ]

    # validate inputs
    if args.select_alpha is not None:
        log.debug(
            "Checking that controls include selected value of alpha %f",
            args.select_alpha,
        )
        if missing_passes := [
            pass_ct
            for pass_ct, controls in enumerate(controls_list, 1)
            if args.select_alpha not in controls.alpha.values
        ]:
            raise ValueError(
                f"Passes {missing_passes} do not summarize selected value of"
                f" alpha={args.select_alpha}"
            )
    if len(controls_list) > 1:
        if args.select_alpha is None:
            log.debug("Checking that controls share same values of alpha")
            controls_alpha = set(controls_list[-1].alpha.values)
            for pass_idx, controls in enumerate(controls_list[:-1]):
                if missing_alpha := controls_alpha.symmetric_difference(
                    controls.alpha.values
                ):
                    raise ValueError(
                        f"Pass {1 + pass_idx} and {len(controls_list)} controls have"
                        f" different alpha (nonoverlapping: {missing_alpha})"
                    )
        log.debug("Checking that controls defined with same experiments")
        controls_experiments = set(controls_list[-1].prefixes)
        for pass_idx, controls in enumerate(controls_list[:-1]):
            if missing := controls_experiments.symmetric_difference(controls.prefixes):
                raise ValueError(
                    f"Pass {1 + pass_idx} and {len(controls_list)} controls have"
                    f" different experiments (nonoverlapping prefixes: {missing}"
                )
    if args.select_alpha_case is not None:
        log.debug(
            "alpha_case is specified: %f",
            args.select_alpha_case,
        )

    # define cases to quantify/compare
    if args.select_cases_prefixes:
        cases_list[-1] = cases_list[-1][args.select_cases_prefixes]
    elif args.drop_cases_prefixes:
        cases_list[-1] = cases_list[-1].drop_prefixes(args.drop_cases_prefixes)
    if len(cases_list) > 1:
        for pass_idx, cases in enumerate(cases_list[:-1]):
            try:
                cases_list[pass_idx] = cases[cases_list[-1].prefixes]
            except KeyError as err:
                missing = set(cases_list[-1].prefixes) - set(cases.prefixes)
                raise ValueError(
                    f"Pass {1 + pass_idx} cases do not have prefixes {missing}"
                    f" from pass {len(cases_list)}"
                ) from err

    if not args.allow_prefix_overlap and (
        overlap := set(controls_list[-1].prefixes) & set(cases_list[-1].prefixes)
    ):
        raise ValueError(f"Controls and cases have overlapping prefixes {overlap}")

    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__
    metadata["n_controls"] = controls_list[-1].num_prefixes
    log.debug("Controls prefixes: %s", controls_list[-1].prefixes)
    metadata["n_cases"] = cases_list[-1].num_prefixes
    log.debug("Cases prefixes: %s", sorted(cases_list[-1].prefixes))
    metadata["abs_dpsi_threshold"] = args.abs_dpsi_threshold
    metadata["select_alpha"] = args.select_alpha  # can be None
    metadata["select_alpha_case"] = args.select_alpha_case  # can be None

    log.debug("Checking that events match and identifying passes with zero events")
    ignore_passes_idx: List[int] = list()
    for pass_idx, (controls, cases) in enumerate(zip(controls_list, cases_list)):
        if not controls.events_df.equals(cases.events_df):
            raise ValueError(
                f"Pass {1 + pass_idx} controls/cases do share identical events"
                f" ({controls.events_df = }, {cases.events_df})"
            )
        if controls.num_events == 0:
            ignore_passes_idx.append(pass_idx)
    # remove passes with zero events
    for pass_idx in reversed(ignore_passes_idx):
        log.info("Ignoring pass %d because it has zero events", 1 + pass_idx)
        del controls_list[pass_idx]
        del cases_list[pass_idx]
        del sg_list[pass_idx]

    log.info(
        "Quantifying potential outliers from cases %s relative to controls %s",
        cases_list,
        controls_list,
    )
    df_list: List[pd.DataFrame] = [
        nm.PsiOutliers(cases, controls, args.select_alpha_case).to_dataframe(
            controls_min_experiments=args.min_experiments,
            alpha=args.select_alpha,
            abs_dpsi_mask_threshold=args.abs_dpsi_threshold,
            show_progress=args.show_progress,
        )
        for cases, controls in zip(cases_list, controls_list)
    ]
    log.info("Annotating quantifications with splicegraph information")
    df: pd.DataFrame
    if len(df_list) > 1:
        # keep any information compatible with current splicegraph, even if it
        # isn't a strict LSV in current splicegraph
        events = sg_list[-1].exon_connections.lsvs(
            nm.constants.SelectLSVs.PERMISSIVE_LSVS
        )
        # get events found in controls
        events_list = [
            controls.get_events(sg.introns, sg.junctions)
            for controls, sg in zip(controls_list, sg_list)
        ]
        # merge df_list onto events
        df = events.merge_dataframes(df_list, events_list, annotated=annotated)
    elif len(df_list) == 1:
        df = (
            controls_list[0]
            .get_events(sg_list[0].introns, sg_list[0].junctions)
            .ec_dataframe(annotated=annotated)
            .join(df_list[0], how="inner", on="ec_idx")
        )
    else:
        df = pd.DataFrame({})

    try:
        output_name = args.output_tsv.name
    except AttributeError:
        output_name = args.output_tsv
    log.info("Saving metadata and table of comparisons to %s", output_name)
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    (
        df
        # manually format probability columns
        .pipe(
            lambda df: df.assign(
                **{
                    col: df[col].apply(lambda x: f"{x:.3e}")
                    for col in df.columns
                    if "probability" in col
                }
            )
        )
        # numeric columns need at most 4 digits precision
        .round(4)
        # save result in TSV format
        .to_csv(args.output_tsv, sep="\t", index=False)
    )

    if args.events_summary or args.genes_summary:
        df_events = nm.PsiOutliers.summarize_df_events(df) if df_list else df
        if args.events_summary:
            try:
                output_name = args.events_summary.name
            except AttributeError:
                output_name = args.events_summary
            log.info(
                "Saving metadata and summary per splicing event to %s", output_name
            )
            args.events_summary.write(
                "# {}\n".format(metadata_json.replace("\n", "\n# "))
            )
            (
                df_events
                # manually format probability columns
                .pipe(
                    lambda df: df.assign(
                        **{
                            col: df[col].apply(lambda x: f"{x:.3e}")
                            for col in df.columns
                            if "probability" in col
                        }
                    )
                )
                # numeric columns need at most 4 digits precision
                .round(4)
                # save result in TSV format
                .to_csv(args.events_summary, sep="\t", index=True)
            )
        if args.genes_summary:
            try:
                output_name = args.genes_summary.name
            except AttributeError:
                output_name = args.genes_summary
            log.info("Saving metadata and summary per gene to %s", output_name)
            args.genes_summary.write(
                "# {}\n".format(metadata_json.replace("\n", "\n# "))
            )
            df_genes = (
                nm.PsiOutliers.summarize_df_genes(df_events) if df_list else df_events
            )
            (
                df_genes
                # manually format probability columns
                .pipe(
                    lambda df: df.assign(
                        **{
                            col: df[col].apply(lambda x: f"{x:.3e}")
                            for col in df.columns
                            if "probability" in col
                        }
                    )
                )
                # numeric columns need at most 4 digits precision
                .round(4)
                # save result in TSV format
                .to_csv(args.genes_summary, sep="\t", index=True)
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
