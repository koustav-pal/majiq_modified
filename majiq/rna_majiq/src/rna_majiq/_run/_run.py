"""
_run.py

Shared types, functions for majiq-tools subcommands

Author: Joseph K Aicher
"""

import argparse
import os
import sys
from hashlib import md5
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, Optional

import dask
from dask.distributed import Client

from rna_majiq import __version__
from rna_majiq.logger import get_logger, setup_logger


class GenericSubcommand(object):
    """Wrap add_args, run with shared setup (especially for logging)"""

    def __init__(
        self,
        DESCRIPTION: str,
        add_args: Callable[[argparse.ArgumentParser], None],
        run: Callable[[argparse.Namespace], None],
    ) -> None:
        self._description = DESCRIPTION
        self._add_args = add_args
        self._run = run
        return

    @property
    def DESCRIPTION(self) -> str:
        return self._description

    def add_args(self, parser: argparse.ArgumentParser) -> None:
        """Add arguments to provided argument parser"""
        parser.add_argument(
            "--version", action="version", version=f"%(prog)s {__version__}"
        )
        self._add_args(parser)
        # add arguments for overwriting existing paths
        parser.add_argument(
            "--overwrite",
            action="store_true",
            default=False,
            help="Ignore/overwrite outputs when output paths already exist",
        )
        # add parameters for logging
        log_args = parser.add_argument_group("logging arguments")
        log_args.add_argument(
            "--logger",
            type=str,
            default=None,
            help="Save logging output to specified file",
        )
        log_args.add_argument(
            "--logfile-only",
            action="store_true",
            default=False,
            help="Only save logging output to file specified by --logger"
            " (default: logging to stderr)."
            " Requires --logger argument to be specified",
        )
        log_args.add_argument(
            "--silent", action="store_true", default=False, help="Silence the logger"
        )
        log_args.add_argument(
            "--debug",
            action="store_true",
            default=False,
            help="Enable detailed logging for debugging purposes",
        )
        return

    def run(self, args: argparse.Namespace) -> None:
        """Run subcommand with parsed arguments"""
        # set up logging
        setup_logger(
            logfile=args.logger,
            silent=args.silent,
            debug=args.debug,
            logfile_only=args.logfile_only,
        )
        log = get_logger()
        # print information about the run
        from rna_majiq._version import version

        log.info(f"new-majiq v{version}")
        command_str = " ".join(sys.argv)
        log.info(f"Command: {command_str}")
        log.info(f"From: {Path(os.getcwd()).resolve()}")
        log.info(
            "\n".join(
                [
                    "Arguments:",
                    "{",
                    *(
                        f" {key} = {value},"
                        for key, value in vars(args).items()
                        if key != "func"
                    ),
                    "}",
                ]
            )
        )
        # set up working directory
        tmp_work_directory: Optional[TemporaryDirectory] = None
        try:
            if args.work_directory:
                log.info(f"Using retained work directory {args.work_directory}")
            else:
                tmp_work_directory = TemporaryDirectory(
                    dir=args.tmp_work_directory_parent,
                    prefix=md5(command_str.encode()).hexdigest(),
                )
                args.work_directory = Path(tmp_work_directory.name).resolve()
                log.info(f"Using temporary work directory {args.work_directory}")
        except AttributeError:
            pass  # work_directory not defined
        # set up dask?
        client: Optional[Client] = None
        if getattr(args, "use_dask", False):
            # set up scheduler
            if args.scheduler_address:
                client = Client(address=args.scheduler_address)
            elif args.scheduler_file:
                client = Client(scheduler_file=args.scheduler_file)
            else:
                dask.config.set(scheduler="threads")
                dask.config.set(num_workers=args.nthreads)
                log.info(f"Using thread pool (nthreads={args.nthreads})")
            if client:
                log.info(f"Using distributed client {client}")
        # add client to args
        args.client = client
        # run subcommand
        try:
            self._run(args)
        except Exception:
            log.exception("Exiting due to exception:")
            sys.exit(-1)
        finally:
            # if dask cluster, clean up after ourselves
            if client:
                client.close()
                del client, args.client
            if tmp_work_directory:
                tmp_work_directory.cleanup()
        # when done running, note that it was successsfully completed!
        log.info("Finished successfully!")
        return
