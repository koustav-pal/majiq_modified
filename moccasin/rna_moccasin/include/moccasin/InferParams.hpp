/**
 * InferParams.hpp
 */
#pragma once

#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

namespace moccasin {

class InferParams {
 public:
  // params = infer_params(uncorrected, factors, passed)
  static constexpr const char name[] = "infer_params";
  static constexpr int ntypes = 2;
  static constexpr int nin = 3;
  static constexpr int nout = 1;
  static constexpr int nargs = nin + nout;
  static constexpr const char signature[] = "(x,n),(f,n),(n)->(x,f)";
  static constexpr const char doc[] = R"mdoc(
  Get coefficients for OLS on uncorrected data given factors

  Parameters
  ----------
  uncorrected: array[float]
    Core dimensions: (experiment, extra)
    Observed values per experiment
  factors: array[float]
    Core dimensions: (experiment, factor)
    Confounding and nonconfounding factors used to model uncorrected data
  passed: array[bool]
    Core dimensions: (experiment)
    Boolean indicator if observation passed and used for modeling

  Returns
  -------
  params: array[float]
    Core dimensions: (extra, factor)
    OLS model parameters such that dot(factors, params, dim="factor") is the
    OLS estimator of uncorrected given passed data
  )mdoc";
  static char types[ntypes * nargs];
  static PyUFuncGenericFunction funcs[ntypes];
};

}  // namespace moccasin
