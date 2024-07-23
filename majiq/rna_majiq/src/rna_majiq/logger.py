"""
logger.py

Wrapper around Python logging for new-majiq

Author: Joseph K Aicher
"""

import logging
from pathlib import Path
from typing import List, Optional, Union


def get_logger(allow_setup: bool = True) -> logging.Logger:
    """Get new-majiq logger

    allow_setup: bool
        Set up logger if it hasn't already been when True
        (set to False when running setup)
    """
    # get logger, which we will add handles to
    logger = logging.getLogger(name="new-majiq")
    if allow_setup and not logger.handlers:
        # if no handlers, then it hasn't been setup yet -- make it silent
        setup_logger(silent=True)
    return logger


def _default_formatter() -> logging.Formatter:
    """How we want log messsages to look"""
    # set up formatting
    log_format = "%(asctime)s (%(levelname)s) - %(message)s"
    formatter = logging.Formatter(log_format)
    return formatter


def setup_logger(
    logfile: Optional[Union[str, Path]] = None,
    silent: bool = False,
    debug: bool = False,
    logfile_only: bool = False,
) -> None:
    """Setup logger to print output to logfile as well as stderr

    Setup logger to print output to logfile as well as stderr. If logfile_only,
    then logfile must be specified, and logging will be written only to logfile
    and not stderr.
    """
    handlers: List[logging.Handler] = [logging.StreamHandler()]
    if logfile_only:
        if not logfile:
            raise ValueError(
                "Cannot request logging to go to logfile only without specifying logfile"
            )
        handlers = []
    if logfile:
        handlers.append(logging.FileHandler(logfile, mode="a"))
    for handler in handlers:
        handler.setFormatter(_default_formatter())
    # get logger, remove existing handlers, add these handlers to it
    logger = get_logger(allow_setup=False)
    for old_handler in logger.handlers:
        logger.removeHandler(old_handler)
    for handler in handlers:
        logger.addHandler(handler)
    # set up level for logging
    if debug:
        logger.setLevel(logging.DEBUG)
    elif silent:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    # done with setup
    return None
