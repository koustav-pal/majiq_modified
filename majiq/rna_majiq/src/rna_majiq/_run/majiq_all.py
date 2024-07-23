"""
majiq_all.py

Single command enumerating all scripts

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from rna_majiq import __version__
from rna_majiq._run._run import GenericSubcommand
from rna_majiq._run.build_combine import subcommand as combine
from rna_majiq._run.build_gff3 import subcommand as gff3
from rna_majiq._run.build_pipeline import subcommand as build_pipeline
from rna_majiq._run.build_update import subcommand as update
from rna_majiq._run.cite import subcommand as cite
from rna_majiq._run.deltapsi import subcommand as deltapsi
from rna_majiq._run.heterogen import subcommand as heterogen
from rna_majiq._run.legacy_majiq import subcommand as legacy_majiq
from rna_majiq._run.legacy_splicegraph import subcommand as legacy_splicegraph
from rna_majiq._run.moccasin import subcommand_coverage_infer as moccasin_coverage_infer
from rna_majiq._run.moccasin import subcommand_coverage_model as moccasin_coverage_model
from rna_majiq._run.moccasin import subcommand_factors_infer as moccasin_factors_infer
from rna_majiq._run.moccasin import subcommand_factors_model as moccasin_factors_model
from rna_majiq._run.moccasin import subcommand_pipeline as moccasin_pipeline
from rna_majiq._run.psi_controls import subcommand as psi_controls
from rna_majiq._run.psi_coverage import subcommand as psi_coverage
from rna_majiq._run.psi_group import subcommand as psi_group
from rna_majiq._run.psi_outliers import subcommand as psi_outliers
from rna_majiq._run.quantify import subcommand as psi
from rna_majiq._run.sg_coverage import subcommand as sg_coverage
from rna_majiq._run.sg_coverage_summarize import subcommand as sg_coverage_summarize
from rna_majiq._run.sj import subcommand as sj

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "build-pipeline": build_pipeline,
    "gff3": gff3,
    "sj": sj,
    "build": update,
    "combine": combine,
    "psi-coverage": psi_coverage,
    "psi-group": psi_group,
    "quantify": psi,
    "deltapsi": deltapsi,
    "heterogen": heterogen,
    "psi-controls": psi_controls,
    "psi-outliers": psi_outliers,
    "moccasin-pipeline": moccasin_pipeline,
    "moccasin-factors-model": moccasin_factors_model,
    "moccasin-factors-infer": moccasin_factors_infer,
    "moccasin-coverage-model": moccasin_coverage_model,
    "moccasin-coverage-infer": moccasin_coverage_infer,
    "sg-coverage": sg_coverage,
    "sg-coverage-summary": sg_coverage_summarize,
    "legacy-majiq": legacy_majiq,
    "legacy-splicegraph": legacy_splicegraph,
    "cite": cite,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="All tools to detect, quantify, and analyze RNA splicing"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
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
