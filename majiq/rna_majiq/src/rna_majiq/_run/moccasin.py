"""
moccasin.py

Perform batch correction using the moccasin algorithm

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict, List, Union

import dask

import rna_majiq as nm
from rna_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    ResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    resources_args,
    select_prefixes_args,
)
from rna_majiq._run._run import GenericSubcommand
from rna_majiq.logger import get_logger


def _args_factors(parser: argparse.ArgumentParser) -> None:
    """add arguments to specify factors to load"""
    factors = parser.add_argument_group("required input confounders arguments")
    factors_ex = factors.add_mutually_exclusive_group(required=True)
    factors_ex.add_argument(
        "--intercept-only",
        dest="intercept_only",
        action="store_true",
        default=None,
        help="Only use non-confounding intercept term"
        " (only makes sense if augmenting with discovered unknown factors)",
    )
    factors_ex.add_argument(
        "--factors-tsv",
        metavar="TSV",
        dest="factors_tsv",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="Path to TSVs with model matrix for all prefixes being processed"
        " (required column: prefix)",
    )
    factors_ex.add_argument(
        "--factors-zarr",
        metavar="ZARR",
        dest="factors_zarr",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="Paths to factors matrices including desired prefixes. If a prefix"
        " is in multiple files, the first listed factors file with that prefix"
        " is used",
    )
    parser.add_argument(
        "--confounding",
        metavar="VAR",
        type=str,
        nargs="+",
        default=list(),
        help="Names of confounding variables, only used with --factors-tsv",
    )
    return


def _get_factors(
    prefix: Union[str, List[str]], args: argparse.Namespace
) -> nm.moccasin.Factors:
    """get factors matrix for specified prefix(es)"""
    factors: nm.moccasin.Factors
    if args.intercept_only:
        factors = nm.moccasin.Factors.intercept_only(prefix)
    elif args.factors_tsv:
        factors = nm.moccasin.Factors.from_tsv(
            args.factors_tsv, args.confounding, prefixes=prefix
        )
    else:
        factors = nm.moccasin.Factors.from_zarr(args.factors_zarr, prefixes=prefix)
    get_logger().info("Using %s\n%s", factors, factors.factor_report())
    return factors


def _args_ruv_parameters(
    parser: argparse.ArgumentParser,
    default_max_new_factors: int = nm.constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
) -> None:
    ruv = parser.add_argument_group("unknown confounders modeling arguments")
    ruv.add_argument(
        "--ruv-max-new-factors",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=default_max_new_factors,
        help="Default number of output new factors for FactorsModel, if any"
        " (default: %(default)s)",
    )
    ruv.add_argument(
        "--ruv-max-events",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        help="Maximum number of events to model residuals for (default: %(default)s)",
    )
    return


def args_factors_model(parser: argparse.ArgumentParser) -> None:
    """arguments for creating model of unknown factors"""
    _args_factors(parser)  # get factors information
    parser.add_argument(
        "factors_model",
        type=ResolvedPath,
        help="Path for output model for unknown confounders",
    )
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths for input psi coverage files for unknown confounders model",
    )
    select_prefixes_args(parser)
    _args_ruv_parameters(parser)
    resources_args(parser, use_dask=True)
    return


def args_factors_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for inferring factors using factors model"""
    _args_factors(parser)
    parser.add_argument(
        "factors_model",
        type=ExistingResolvedPath,
        help="Path for input model of unknown confounders",
    )
    parser.add_argument(
        "output",
        type=ResolvedPath,
        help="Path for output factors (known and unknown)",
    )
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to psi coverage files for which known/unknown factors will be saved",
    )
    parser.add_argument(
        "--num-unknown-factors",
        type=check_nonnegative_factory(int, reject_zero=False),
        default=None,
        help="If specified, override number of output unknown confounders"
        " produced by model (default: set in majiq-moccasin factors-model)",
    )
    select_prefixes_args(parser)
    resources_args(parser, use_dask=True)
    return


def args_coverage_model(parser: argparse.ArgumentParser) -> None:
    """arguments for modeling coverage"""
    parser.add_argument(
        "coverage_model",
        type=ResolvedPath,
        help="Output path for coverage model",
    )
    _args_factors(parser)
    parser.add_argument(
        "psicov",
        nargs="+",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        help="Paths for input psi coverage files for coverage model",
    )
    parser.add_argument(
        "--prefix-nchunks",
        type=int,
        default=1,
        help="Maximum number of prefixes loaded per chunk of PsiCoverage."
        " This speeds up reduction over prefixes performed when calculating"
        " gramian and projection matrices to solve model parameters."
        " This coarser chunking can increase memory usage and is only done"
        " within individual input PsiCoverage files (and not between them)."
        " Non-positive values load all prefixes per file together."
        " (default: %(default)s)",
    )
    select_prefixes_args(parser)
    resources_args(parser, use_dask=True)
    return


def args_coverage_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for getting corrected lsv coverage"""
    _args_factors(parser)
    parser.add_argument(
        "coverage_model",
        type=ExistingResolvedPath,
        help="Input path for coverage model",
    )
    parser.add_argument(
        "corrected_psicov",
        type=ResolvedPath,
        help="Path for output corrected psi coverage",
    )
    parser.add_argument(
        "original_psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to original, uncorrected psi coverage",
    )
    select_prefixes_args(parser)
    resources_args(parser, use_dask=True)
    return


def run_factors_model(args: argparse.Namespace) -> None:
    """model unknown factors using input coverage"""
    if not args.overwrite and args.factors_model.exists():
        raise ValueError(
            f"Output factors model {args.factors_model} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Opening coverage from %d PSI coverage files", len(args.psicov))
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psicov = psicov[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psicov = psicov.drop_prefixes(args.drop_prefixes)
    log.info("Setting up model matrix of known factors")
    factors: nm.moccasin.Factors = _get_factors(psicov.prefixes, args)
    log.info("Learning model for unknown confounding factors from %s", psicov)
    log.debug("Prefixes: %s", psicov.prefixes)
    model = nm.moccasin.FactorsModel.train(
        psicov,
        factors,
        default_num_factors=args.ruv_max_new_factors,
        max_events=args.ruv_max_events,
        log=log,
        show_progress=args.show_progress,
    )
    # report about explained variance
    if model.max_num_factors > 0:
        log.info(model.explained_variance_report())
    # save model
    log.info("Saving model to %s", args.factors_model)
    model.to_zarr(args.factors_model)
    return


def run_factors_infer(args: argparse.Namespace):
    """compute unknown factors using input coverage"""
    if not args.overwrite and args.output.exists():
        raise ValueError(
            f"Output factors {args.output} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Opening coverage from %d PsiCoverage files", len(args.psicov))
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psicov = psicov[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psicov = psicov.drop_prefixes(args.drop_prefixes)
    log.info("Inferring factors for %s", psicov)
    log.debug("Prefixes: %s", psicov.prefixes)
    factors = _get_factors(psicov.prefixes, args)
    log.info("Using known factors %s", factors)
    log.info("Loading model for unknown factors from %s", args.factors_model)
    model = nm.moccasin.FactorsModel.from_zarr(args.factors_model)
    log.info("Solving for unknown confounders")
    updated_factors = model.predict(
        psicov,
        factors,
        num_unknown_factors=args.num_unknown_factors,
        show_progress=args.show_progress,
    )
    log.info("Saving %s to %s", updated_factors, args.output)
    updated_factors.to_zarr(args.output)
    return


def run_coverage_model(args: argparse.Namespace) -> None:
    """Determine parameters for coverage model given factors"""
    if not args.overwrite and args.coverage_model.exists():
        raise ValueError(
            f"Output coverage model {args.coverage_model} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Opening coverage from %d PsiCoverage files", len(args.psicov))
    psicov = nm.PsiCoverage.from_zarr(
        args.psicov, ec_idx_nchunks=1, prefix_nchunks=args.prefix_nchunks
    )
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psicov = psicov[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psicov = psicov.drop_prefixes(args.drop_prefixes)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args)
    log.info("Solving for model parameters from %s using %s", psicov, factors)
    log.debug("Prefixes: %s", psicov.prefixes)
    model = nm.moccasin.CoverageModel.train(psicov, factors)
    log.info("Saving model parameters to %s", args.coverage_model)
    with dask.config.set(
        {"optimization.fuse.ave-width": 2 * len(psicov.df.chunks["ec_idx"])}
    ):
        model.to_zarr(args.coverage_model, show_progress=args.show_progress)
    return


def run_coverage_infer(args: argparse.Namespace) -> None:
    """Create corrected LSV coverage file"""
    if not args.overwrite and args.corrected_psicov.exists():
        raise ValueError(
            f"Output corrected coverage {args.corrected_psicov} already exists"
            " (add `--overwrite` to force replacing it)"
        )
    log = get_logger()
    log.info("Opening coverage from %d PsiCoverage files", len(args.original_psicov))
    psicov = nm.PsiCoverage.from_zarr(
        args.original_psicov, ec_idx_nchunks=1, prefix_nchunks=1
    )
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psicov = psicov[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psicov = psicov.drop_prefixes(args.drop_prefixes)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args)
    log.info("Opening up model parameters from %s", args.coverage_model)
    model = nm.moccasin.CoverageModel.from_zarr(args.coverage_model)
    log.info(
        "Updating %s using %s and model from %s", psicov, factors, args.coverage_model
    )
    log.debug("Prefixes: %s", psicov.prefixes)
    psicov = model.predict(psicov, factors, model_path=str(args.coverage_model))
    log.info("Saving corrected coverage to %s", args.corrected_psicov)
    with dask.config.set(
        {"optimization.fuse.ave-width": 2 * len(psicov.df.chunks["ec_idx"])}
    ):
        psicov.to_zarr(args.corrected_psicov, show_progress=args.show_progress)
    return


def args_pipeline(parser: argparse.ArgumentParser) -> None:
    """arguments for pipeline through model/inference steps of factors/coverage"""
    _args_factors(parser)
    _args_ruv_parameters(parser, default_max_new_factors=0)
    parser.add_argument(
        "output_dir",
        type=NewResolvedPath,
        help="Path for new directory for output files",
    )
    parser.add_argument(
        "psicov",
        nargs="+",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        help="Paths for input psi coverage files to perform batch correction on",
    )
    select_prefixes_args(parser)
    factor_prefixes = parser.add_argument_group("factor modeling experiments arguments")
    factor_prefixes_ex = factor_prefixes.add_mutually_exclusive_group()
    factor_prefixes_ex.add_argument(
        "--factors-prefixes-include",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly include in factors-model step"
        " (default: include all)",
    )
    factor_prefixes_ex.add_argument(
        "--factors-prefixes-exclude",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly exclude in factors-model step,"
        " i.e. use all others (default: include all)",
    )
    coverage_prefixes = parser.add_argument_group(
        "coverage modeling experiments arguments"
    )
    coverage_prefixes_ex = coverage_prefixes.add_mutually_exclusive_group()
    coverage_prefixes_ex.add_argument(
        "--coverage-prefixes-include",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly include in coverage-model step"
        " (default: include all)",
    )
    coverage_prefixes_ex.add_argument(
        "--coverage-prefixes-exclude",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly exclude in coverage-model step,"
        " i.e. use all others (default: include all)",
    )
    resources_args(parser, use_dask=True)
    return


def run_pipeline(args: argparse.Namespace) -> None:
    """model factors, get updated factors, model coverage, get updated coverage"""
    log = get_logger()
    log.info("Opening coverage from %d PsiCoverage files", len(args.psicov))
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    if args.select_prefixes:
        log.debug("Selecting specified prefixes")
        psicov = psicov[args.select_prefixes]
    elif args.drop_prefixes:
        log.debug("Dropping specified prefixes")
        psicov = psicov.drop_prefixes(args.drop_prefixes)
    psicov_subset: nm.PsiCoverage  # used for modeling steps
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args)
    log.info("Modeling confounders pipeline for %s starting with %s", psicov, factors)

    if args.ruv_max_new_factors > 0:
        if args.factors_prefixes_include:
            log.info("Modeling factors using specified prefixes")
            psicov_subset = psicov[args.factors_prefixes_include]
        elif args.factors_prefixes_exclude:
            log.info("Modeling factors excluding specified prefixes")
            psicov_subset = psicov.drop_prefixes(args.factors_prefixes_exclude)
        else:
            psicov_subset = psicov
        log.info(
            "Learning model for unknown confounding factors from %s", psicov_subset
        )
        factors_model = nm.moccasin.FactorsModel.train(
            psicov_subset,
            factors,
            default_num_factors=args.ruv_max_new_factors,
            max_events=args.ruv_max_events,
            log=log,
            show_progress=args.show_progress,
        )
        # report about explained variance
        if factors_model.max_num_factors > 0:
            log.info(factors_model.explained_variance_report())
        # save factors model
        factors_model_path = args.output_dir / "factors_model.zarr"
        log.info("Saving factors model to %s", factors_model_path)
        factors_model.to_zarr(factors_model_path)

        # infer factors for all experiments
        log.info("Updating factors for %s", psicov)
        factors = factors_model.predict(
            psicov, factors, show_progress=args.show_progress
        )
        factors_path = args.output_dir / "factors.zarr"
        log.info("Saving combined known/unknown factors to %s", factors_path)
        factors.to_zarr(factors_path)
    # done modeling/updating factors if modeling unknown factors

    # for coverage modeling/inference steps, use chunks on ec_idx and prefixes
    psicov = nm.PsiCoverage.from_zarr(args.psicov, ec_idx_nchunks=1, prefix_nchunks=1)
    if args.coverage_prefixes_include:
        log.info("Modeling coverage using specified prefixes")
        psicov_subset = psicov[args.coverage_prefixes_include]
    elif args.coverage_prefixes_exclude:
        log.info("Modeling coverage excluding specified prefixes")
        psicov_subset = psicov.drop_prefixes(args.coverage_prefixes_exclude)
    else:
        psicov_subset = psicov
    log.info("Learning model for observed coverage from %s", psicov_subset)
    coverage_model = nm.moccasin.CoverageModel.train(psicov_subset, factors)
    models_path = args.output_dir / "coverage_model.zarr"
    log.info("Saving coverage model to %s", models_path)
    coverage_model.to_zarr(models_path, show_progress=args.show_progress)
    # reload coverage_model from disk
    coverage_model = nm.moccasin.CoverageModel.from_zarr(models_path)

    log.info("Updating coverage for %s", psicov)
    psicov = nm.PsiCoverage.from_zarr(args.psicov, ec_idx_nchunks=1, prefix_nchunks=1)
    psicov = coverage_model.predict(psicov, factors, model_path=str(models_path))
    psicov_path = args.output_dir / "corrected.psicov"
    log.info("Saving corrected coverage to %s", psicov_path)
    psicov.to_zarr(psicov_path, show_progress=args.show_progress)
    return


subcommand_factors_model = GenericSubcommand(
    "Build model of unknown confounding factors using input PsiCoverage",
    args_factors_model,
    run_factors_model,
)
subcommand_factors_infer = GenericSubcommand(
    "Use unknown confounding factors model on inputs to infer known/unknown factors",
    args_factors_infer,
    run_factors_infer,
)
subcommand_coverage_model = GenericSubcommand(
    "Build model of PsiCoverage using input coverage and factors",
    args_coverage_model,
    run_coverage_model,
)
subcommand_coverage_infer = GenericSubcommand(
    "Use PsiCoverage model to infer corrected coverage given inputs",
    args_coverage_infer,
    run_coverage_infer,
)
subcommand_pipeline = GenericSubcommand(
    "majiq-moccasin pipeline for PsiCoverage batch correction with"
    " known/unknown factors",
    args_pipeline,
    run_pipeline,
)


SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "factors-model": subcommand_factors_model,
    "factors-infer": subcommand_factors_infer,
    "coverage-model": subcommand_coverage_model,
    "coverage-infer": subcommand_coverage_infer,
    "pipeline": subcommand_pipeline,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to model factors and PsiCoverage using MOCCASIN"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {nm.__version__}"
    )
    # add subparsers
    subparsers = parser.add_subparsers(required=True, help="")
    for src_name, src_module in SUBPARSER_SOURCES.items():
        src_parser = subparsers.add_parser(
            src_name,
            help=src_module.DESCRIPTION,
            description=src_module.DESCRIPTION,
        )
        src_parser.set_defaults(func=src_module.run)
        src_module.add_args(src_parser)

    # check length of input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # parse arguments now
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
