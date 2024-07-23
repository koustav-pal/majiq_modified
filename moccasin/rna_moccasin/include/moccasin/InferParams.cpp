/**
 * InferParams.cpp
 */

#include "InferParams.hpp"

#include <Eigen/QR>
#include <limits>
#include <vector>

#include "NumpyEigen.hpp"
#include "NumpyMaskToIndices.hpp"

namespace moccasin {

template <typename RealT, NumpyStrideType t0 = NumpyStrideType::UNDETERMINED,
          NumpyStrideType t1 = NumpyStrideType::UNDETERMINED,
          NumpyStrideType t2 = NumpyStrideType::UNDETERMINED,
          NumpyStrideType t3 = NumpyStrideType::UNDETERMINED>
void InferParamsOuterDispatch(char** args, const npy_intp* dims,
                              const npy_intp* strides, void* data);

char InferParams::types[] = {
    // npy_float function
    NPY_FLOAT,
    NPY_FLOAT,
    NPY_BOOL,
    NPY_FLOAT,
    // npy_double function
    NPY_DOUBLE,
    NPY_DOUBLE,
    NPY_BOOL,
    NPY_DOUBLE,
};
PyUFuncGenericFunction InferParams::funcs[] = {
    InferParamsOuterDispatch<npy_float>,
    InferParamsOuterDispatch<npy_double>,
};

template <typename RealT, NumpyStrideType t0, NumpyStrideType t1,
          NumpyStrideType t2, NumpyStrideType t3>
void InferParamsOuterImpl(char** args, const npy_intp* dims,
                          const npy_intp* strides, void* data) {
  const npy_intp* core_dims = 1 + dims;
  const npy_intp* core_strides = InferParams::nargs + strides;

  auto dim_x = core_dims[0];
  auto dim_n = core_dims[1];
  auto dim_f = core_dims[2];

  auto s_uncorrected_x = core_strides[0];
  auto s_uncorrected_n = core_strides[1];
  if constexpr (t0 == NumpyStrideType::CORE_ALIGNED) {
    s_uncorrected_x /= sizeof(RealT);
    s_uncorrected_n /= sizeof(RealT);
  }
  auto s_factors_f = core_strides[2];
  auto s_factors_n = core_strides[3];
  if constexpr (t1 == NumpyStrideType::CORE_ALIGNED) {
    s_factors_f /= sizeof(RealT);
    s_factors_n /= sizeof(RealT);
  }
  auto s_passed_n = core_strides[4];
  auto s_params_x = core_strides[5];
  auto s_params_f = core_strides[6];
  if constexpr (t3 == NumpyStrideType::CORE_ALIGNED) {
    s_params_x /= sizeof(RealT);
    s_params_f /= sizeof(RealT);
  }

  auto uncorrected = InitNumpyArray<RealT, t0>(args[0], dim_x, s_uncorrected_x,
                                               dim_n, s_uncorrected_n);
  auto factors = InitNumpyArray<RealT, t1>(args[1], dim_f, s_factors_f, dim_n,
                                           s_factors_n);
  std::vector<int> passed_indices;
  passed_indices.reserve(dim_n);
  bool passed_insufficient = false;
  if constexpr (t2 == NumpyStrideType::BROADCAST_ZERO) {
    UpdateIndices(passed_indices, args[2], dim_n, s_passed_n);
    passed_insufficient = static_cast<npy_intp>(passed_indices.size()) < dim_f;
  }
  auto params =
      InitNumpyArray<RealT, t3>(args[3], dim_x, s_params_x, dim_f, s_params_f);

  constexpr bool uncorrected_const = t0 == NumpyStrideType::BROADCAST_ZERO;
  constexpr bool factors_const = t1 == NumpyStrideType::BROADCAST_ZERO;
  constexpr bool passed_const = t2 == NumpyStrideType::BROADCAST_ZERO;

  using QR_t = Eigen::ColPivHouseholderQR<MatrixT<RealT>>;
  QR_t qr{};
  constexpr bool QR_const = factors_const && passed_const;
  if constexpr (QR_const) {
    if (!passed_insufficient) {
      qr.compute(factors(passed_indices, Eigen::all));
    }
  }

  for (npy_intp i_broadcast = 0; i_broadcast < *dims; ++i_broadcast) {
    // update args
    if constexpr (t0 != NumpyStrideType::BROADCAST_ZERO) {
      ReadNumpyArray(uncorrected, args[0], dim_x, s_uncorrected_x, dim_n,
                     s_uncorrected_n);
    }
    if constexpr (t1 != NumpyStrideType::BROADCAST_ZERO) {
      ReadNumpyArray(factors, args[1], dim_f, s_factors_f, dim_n, s_factors_n);
    }
    if constexpr (t2 != NumpyStrideType::BROADCAST_ZERO) {
      UpdateIndices(passed_indices, args[2], dim_n, s_passed_n);
      passed_insufficient =
          static_cast<npy_intp>(passed_indices.size()) < dim_f;
    }
    if constexpr (t3 == NumpyStrideType::CORE_ALIGNED) {
      // we need to update the map location, but otherwise we don't need to
      // read the output values
      ReadNumpyArray(params, args[3], dim_x, s_params_x, dim_f, s_params_f);
    }
    if (!passed_insufficient) {
      // factors(n, f) * params(f, x) = uncorrected(n, x)
      // QR decomposition of factors selecting passed indices
      if constexpr (!QR_const) {
        qr.compute(factors(passed_indices, Eigen::all));
      }
      if (qr.rank() == dim_f) {
        params = qr.solve(uncorrected(passed_indices, Eigen::all));
      } else {
        params.fill(std::numeric_limits<RealT>::quiet_NaN());
      }
    } else {
      params.fill(std::numeric_limits<RealT>::quiet_NaN());
    }
    // write params back to numpy buffer
    WriteNumpyArray(params, args[3], dim_x, s_params_x, dim_f, s_params_f);
    // update args with next strides, if necessary
    if constexpr (t0 != NumpyStrideType::BROADCAST_ZERO) {
      args[0] += strides[0];
    }
    if constexpr (t1 != NumpyStrideType::BROADCAST_ZERO) {
      args[1] += strides[1];
    }
    if constexpr (t2 != NumpyStrideType::BROADCAST_ZERO) {
      args[2] += strides[2];
    }
    if constexpr (t3 != NumpyStrideType::BROADCAST_ZERO) {
      args[3] += strides[3];
    }
  }
}

template <typename RealT, NumpyStrideType t0, NumpyStrideType t1,
          NumpyStrideType t2, NumpyStrideType t3>
void InferParamsOuterDispatch(char** args, const npy_intp* dims,
                              const npy_intp* strides, void* data) {
  if constexpr (t0 == NumpyStrideType::UNDETERMINED) {
    switch (GetNumpyStrideType<RealT>(strides[0], strides[4], strides[5])) {
      case NumpyStrideType::BROADCAST_ZERO:
        return InferParamsOuterDispatch<RealT, NumpyStrideType::BROADCAST_ZERO,
                                        t1, t2, t3>(args, dims, strides, data);
      case NumpyStrideType::CORE_ALIGNED:
        return InferParamsOuterDispatch<RealT, NumpyStrideType::CORE_ALIGNED,
                                        t1, t2, t3>(args, dims, strides, data);
      case NumpyStrideType::CORE_UNALIGNED:
        return InferParamsOuterDispatch<RealT, NumpyStrideType::CORE_UNALIGNED,
                                        t1, t2, t3>(args, dims, strides, data);
      default:
        throw std::runtime_error(
            "GetNumpyStrideType returned unexpected value");
    }
  }
  if constexpr (t1 == NumpyStrideType::UNDETERMINED) {
    switch (GetNumpyStrideType<RealT>(strides[1], strides[6], strides[7])) {
      case NumpyStrideType::BROADCAST_ZERO:
        return InferParamsOuterDispatch<
            RealT, t0, NumpyStrideType::BROADCAST_ZERO, t2, t3>(args, dims,
                                                                strides, data);
      case NumpyStrideType::CORE_ALIGNED:
        return InferParamsOuterDispatch<RealT, t0,
                                        NumpyStrideType::CORE_ALIGNED, t2, t3>(
            args, dims, strides, data);
      case NumpyStrideType::CORE_UNALIGNED:
        return InferParamsOuterDispatch<
            RealT, t0, NumpyStrideType::CORE_UNALIGNED, t2, t3>(args, dims,
                                                                strides, data);
      default:
        throw std::runtime_error(
            "GetNumpyStrideType returned unexpected value");
    }
  }
  if constexpr (t2 == NumpyStrideType::UNDETERMINED) {
    switch (GetNumpyStrideType<bool>(strides[2], strides[8])) {
      case NumpyStrideType::BROADCAST_ZERO:
        return InferParamsOuterDispatch<RealT, t0, t1,
                                        NumpyStrideType::BROADCAST_ZERO, t3>(
            args, dims, strides, data);
      case NumpyStrideType::CORE_ALIGNED:
        return InferParamsOuterDispatch<RealT, t0, t1,
                                        NumpyStrideType::CORE_ALIGNED, t3>(
            args, dims, strides, data);
      case NumpyStrideType::CORE_UNALIGNED:
        return InferParamsOuterDispatch<RealT, t0, t1,
                                        NumpyStrideType::CORE_UNALIGNED, t3>(
            args, dims, strides, data);
      default:
        throw std::runtime_error(
            "GetNumpyStrideType returned unexpected value");
    }
  }
  if constexpr (t3 == NumpyStrideType::UNDETERMINED) {
    switch (GetNumpyStrideType<RealT>(strides[3], strides[9], strides[10])) {
      case NumpyStrideType::CORE_ALIGNED:
        return InferParamsOuterDispatch<RealT, t0, t1, t2,
                                        NumpyStrideType::CORE_ALIGNED>(
            args, dims, strides, data);
      case NumpyStrideType::UNDETERMINED:
        throw std::runtime_error(
            "GetNumpyStrideType returned unexpected value");
      default:
        return InferParamsOuterDispatch<RealT, t0, t1, t2,
                                        NumpyStrideType::CORE_UNALIGNED>(
            args, dims, strides, data);
    }
  }
  return InferParamsOuterImpl<RealT, t0, t1, t2, t3>(args, dims, strides, data);
}

}  // namespace moccasin
