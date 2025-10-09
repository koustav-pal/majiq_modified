.. _moccasin:

############################
Correcting for batch effects
############################

.. role:: bash(code)
   :language: bash

Correcting for confounding variables/batch effects may be tricky...

.. image:: https://imgs.xkcd.com/comics/confounding_variables.png
   :alt: xkcd on confounding variables
   :target: https://xkcd.com/2560/

... but we should still try.

MOCCASIN is a framework for splicing confounder adjustment [SlaffRadens2021]_.
MAJIQ provides commands (:bash:`majiq moccasin`, subcommands from
:bash:`majiq-moccasin`) to use this framework to correct batch effects in
PsiCoverage files.
The MAJIQ MOCCASIN pipeline can be run on a single machine using the command
:bash:`majiq moccasin`.
This command chains together the following steps:

1. (Optional) model/infer unknown confounding factors from known factors
   (:bash:`majiq-moccasin factors-model`, :bash:`majiq-moccasin factors-infer`)
2. Model raw and bootstrapped coverage per LSV given all factors
   (:bash:`majiq-moccasin coverage-model`)
3. Get corrected raw and bootstrapped coverage per LSV, removing the effect of
   confounding factors (:bash:`majiq-moccasin coverage-infer`).


The MOCCASIN pipeline
=====================

From the perspective of input and output files, MOCCASIN takes as input and
produces as output the following:

Input files
-----------

- Uncorrected PsiCoverage from the MAJIQ builder.
- Observed confounding nonconfounding factors per experiment used to model
  unobserved confounding factors and/or coverage.

Output files
------------

- Model of unobserved confounding factors (`factors_model.zarr`).
- Combined input observed factors and inferred unobserved confounding factors
  for input experiments (`factors.zarr`).
- Model of raw and bootstrapped coverage per LSV (`coverage_model.zarr`).
- Corrected PsiCoverage removing the effect of confounding factors
  (`corrected.psicov`).


Command-line interface
----------------------

From the persepctive of the command-line, the MOCCASIN pipeline specifies the
locations of these inputs and outputs using the required positional arguments:

- :bash:`output_dir`: path for new output directory to which models and
  inferred factors/coverage are saved.
- :bash:`psicov [psicov ...]`: paths for input PsiCoverage files for which
  batch correction is performed.

MOCCASIN also requires input confunding factors to be specified with one
of three mutually-exclusive arguments:

- :bash:`--factors-tsv TSV`: TSV with model matrix of observed factors.
  Required column `prefix` indicates experiment names in input PsiCoverage
  files per row; additional columns are observed factors.
  The argument :bash:`--confounding VAR [VAR ...]` indicates which columns
  correspond to confounding factors (others are treated as nonconfounding).

    - **Note**: MOCCASIN does not automatically add a nonconfounding intercept
      term. You probably want this.
    - **Note**: factors are treated as numeric values.
      Categorical factors need to be input as dummy variables.
      Picking appropriate contrasts for confounding variables is important
      because batch correction involves predictions with the confounding
      variables set to zero.

- :bash:`--factors-zarr ZARR [ZARR ...]`: Paths to factors matrices including
  prefixes from input PsiCoverage flies.
  This is typically output from :bash:`majiq-moccasin` itself.
- :bash:`--intercept-only`: The only observed factor should be a
  non-confounding intercept term.
  This is commonly used with unobserved confounding factors.

As an example, adding the options
:bash:`--factors-tsv {TSV} --confounding from_lane2 from_lane3`
with the following table models a nonconfounding intercept term and corrects
experiments from lane 2 and lane 3:

====== ========= ========== ==========
prefix intercept from_lane2 from_lane3
====== ========= ========== ==========
Xa     1         0          0
Xb     1         1          0
Xc     1         0          1
Xd     1         0          1
...    ...       ...        ...
Ya     1         0          0
Yb     1         0          1
Za     1         1          0
Zb     1         0          1
Zc     1         1          0
====== ========= ========== ==========

To enable modeling/inference of unobserved confounding factors, set
:bash:`--ruv-max-new-factors`.
More detailed explanations of these parameters (and others) can be found by
running :bash:`majiq moccasin --help`.


Tips for setting up your model matrix
--------------------------------------
MOCCASIN learns a linear model per junction or intron and calculates
corrections by setting the confounding factors to zero. (For full details,
please see the Moccasin paper: [SlaffRadens2021]_.) This leads to some considerations for
forming the model matrix:

- Confounder effects are modeled linearly. Consider transforming factor values
  if you believe some re-scaling of the factor is better for linear modeling
  (e.g. log or exponential).
- Often, it makes sense to break a factor into discrete categories, as in the
  lane 1, 2,or 3 example shown above. Then the factor values can be 1 (yes)
  or 0 (no) in each column.
- MOCCASIN requires the model matrix to be full-rank. A consequence is that if you
  have discrete categories, you cannot include all the category columns together with
  the intercept, since then the sum of the category columns would be all 1s, equal to
  the intercept (hence, not full-rank). That is why the lane 1, 2, 3 example
  above leaves out lane 1. Moreover, the left-out factor from such a group is, in effect,
  modeled by the intercept. Since MOCCASIN calculates the corrected values by setting
  all confounding factors to zero, the intercept (hence, the left-out factor) is
  the one represented in the corrected data. So, in the lane 1, 2, or 3 example above,
  the data is corrected to be like all samples were in lane 1.
- Since MOCCASIN calculates corrected values by setting all confounding factors to zero,
  a zero value for each factor should represent the case you want to see in the
  corrected data. For example, suppose you want to include RNA Integrity Number (RIN)
  as a confounding factor. RIN goes from 1 to 10, with 10 the highest-integrity RNA.
  If not transformed, MOCCASIN will model RIN and then correct all samples to RIN = 0,
  which represents off-the-scale degraded RNA. On the other hand, if you replace RIN with
  10-RIN, 0 represents the highest-integrity RNA. Another option is to bucket RIN into
  discrete groups, e.g. high, middle, and low RNA integrity (1/0 values). Leaving the
  "high" category out of the model matrix would lead to MOCCASIN correcting the data be like
  high-integrity RNA.



Individual steps with :bash:`majiq-moccasin`
============================================

The MOCCASIN pipeline chains together 4 different commands from
:bash:`majiq-moccasin`.
There are many cases where you might want to use these commands directly rather
than the pipeline.

These include but are not limited to:

- Creating factors files faster (:bash:`majiq-moccasin factors-infer`): the
  pipeline infers factors for all input experiments on a single machine. With
  large numbers of experiments, it can sometimes be faster to batch the
  experiments to produce multiple output factors files (:bash:`--factors-zarr`
  permits multiple factors files to permit this batching)
- Grouping PsiCoverage files/creating them faster
  (:bash:`majiq-moccasin coverage-infer`):
  the pipeline infers coverage for all input experiments on a single machine to
  a single output file.
  If you wish to maintain groupings from the builder, you need to run
  :bash:`majiq-moccasin coverage-infer` separately on each PsiCoverage file.
  This also allows this step to be run in a distributed fashion.
- Inferring factors/batch correction for new additional experiments
  (i.e. for :ref:`two-pass build <quick-twopass>`).
  While :bash:`majiq moccasin` has flags to include/exclude specific prefixes
  from modeling steps, it still requires the PsiCoverage files to be present at
  the time of modeling.
  This allows for the MOCCASIN models to be applied to new data for additional
  experiments.
