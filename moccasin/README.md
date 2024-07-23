# MOCCASIN

MOCCASIN is a Python package for RNA splicing confounding factor adjustment
from short-read RNA-seq data.
MOCCASIN was designed for use with [MAJIQ][majiq-webpage] but has an API that
should let it work more generally with other kinds of count data.
After correcting for read depth, MOCCASIN models the variation in junction read
rates due to confounding factors and covariates.
MOCCASIN then adjusts the read rates to remove variation due to confounding
factors and outputs updated tables of junction read rates.
Please see [Slaff et al. 2021][moccasin-paper] for details.

   [majiq-webpage]: (https://majiq.biociphers.org/)
   [moccasin-paper]: (https://doi.org/10.1038/s41467-021-23608-9)


This branch is a rewrite of the core functionality of MOCCASIN for increased
efficiency.
At this time, there is no command-line interface on this branch.
Once there is, we will replace this paragraph with a concise explanation of
how MOCCASIN is run on data from MAJIQ (explained in more details in further
sections of this README).

MOCCASIN is open source under the MIT license (see `LICENSE.md`).

Please contact yosephb@upenn.edu with comments or questions.


## Installation

Install using pip (`pip install .`)


### For development

We use black, flake8, and mypy to enforce a uniform code-style and check for
typing.
This is done through [pre-commit hooks][precommit].
This requires `pre-commit` to be set up.
Either use the conda environment file (`environment.yml`) or install using pip
(`pip install pre-commit`).
Then, run `pre-commit install` to integrate it as a hook whenever attempting a
commit to git.
This only checks changes being committed.
It is occasionally worth checking all files with these checks -- this can be
done by running `pre-commit run --all-files`.

   [precommit]: https://pre-commit.com/


## Running MOCCASIN with coverage from MAJIQ v2.2

TODO...
Input files: majiq files, model matrix...
What does the model matrix look like...
Reserved/required column names for model matrix specification
Expected output files...
How to specify output files...
Other required arguments...
How to configure discovery of unknown confounders...
Tell people to run `moccasin --help` for more details when appropriate but
please don't copy the text of the help message here because otherwise we need
to keep it up to date whenever we update that...


## MOCCASIN API

TODO...
Could probably set up sphinx if we set up comments right...


## Technical notes

Some comments about things that might seem weird?
Such as updating duplicate coverage rows independently depending on whether
they belong to LSV
