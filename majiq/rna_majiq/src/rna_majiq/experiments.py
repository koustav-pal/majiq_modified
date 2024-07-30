"""
experiments.py

We infer experiment name information from input BAM names

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Union

import numpy as np
import numpy.typing as npt
import xarray as xr


def bam_experiment_name(p: Union[str, Path]) -> str:
    # return file name with extension removed
    return Path(p).name.rsplit(".", 1)[0]


def min_experiments(
    min_experiments_f: Union[xr.DataArray, npt.ArrayLike],
    num_experiments: int,
):
    if isinstance(min_experiments_f, xr.DataArray):
        min_experiments_f = min_experiments_f.where(
            min_experiments_f >= 1, min_experiments_f * num_experiments
        )
    else:
        min_experiments_f = np.asarray(min_experiments_f)
        min_experiments_f = np.where(
            min_experiments_f >= 1,
            min_experiments_f,
            min_experiments_f * num_experiments,
        )
    return np.clip(min_experiments_f, 1, num_experiments)
