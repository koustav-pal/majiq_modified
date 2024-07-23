/**
 * NumpyEigen.hpp
 *
 * View or copy NumPy arrays as provided by the (g)ufunc interface as Eigen
 * objects (transposed as Eigen is column order by default
 */
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <type_traits>

namespace moccasin {

enum class NumpyStrideType : uint8_t {
  UNDETERMINED,
  /// broadcast stride is zero, so the array is constant
  BROADCAST_ZERO,
  /// core stride(s) are a multiple of the underlying type (can use Eigen::Map)
  CORE_ALIGNED,
  /// core stride(s) are not aligned, so need to copy
  CORE_UNALIGNED
};

template <typename RealT, typename IndexT, typename... OtherT>
inline constexpr bool IsStridesUnaligned(IndexT stride,
                                         OtherT... other_stride) {
  bool stride_unaligned = stride % sizeof(RealT);
  if constexpr (sizeof...(OtherT)) {
    return stride_unaligned || IsStridesUnaligned<RealT>(other_stride...);
  } else {
    return stride_unaligned;
  }
}

template <typename RealT, typename IndexT, typename... OtherT>
inline constexpr NumpyStrideType GetNumpyStrideType(IndexT broadcast_stride,
                                                    IndexT stride,
                                                    OtherT... other_stride) {
  if (broadcast_stride) {
    if (IsStridesUnaligned<RealT>(stride, other_stride...)) {
      return NumpyStrideType::CORE_UNALIGNED;
    } else {
      return NumpyStrideType::CORE_ALIGNED;
    }
  } else {
    return NumpyStrideType::BROADCAST_ZERO;
  }
}

using MatrixStrideT = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
using VectorStrideT = Eigen::InnerStride<Eigen::Dynamic>;

template <typename RealT>
using MatrixT =
    Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <typename RealT>
using VectorT = Eigen::Matrix<RealT, Eigen::Dynamic, 1, Eigen::ColMajor>;

template <typename RealT>
using MatrixMapT = Eigen::Map<MatrixT<RealT>, Eigen::Unaligned, MatrixStrideT>;
template <typename RealT>
using VectorMapT = Eigen::Map<VectorT<RealT>, Eigen::Unaligned, VectorStrideT>;

template <typename RealT, NumpyStrideType stride_type>
using NumpyMatrixT =
    std::conditional_t<stride_type == NumpyStrideType::CORE_ALIGNED,
                       MatrixMapT<RealT>, MatrixT<RealT>>;
template <typename RealT, NumpyStrideType stride_type>
using NumpyVectorT =
    std::conditional_t<stride_type == NumpyStrideType::CORE_ALIGNED,
                       VectorMapT<RealT>, VectorT<RealT>>;

/// Updates m to represent the transpose of the numpy array
template <typename RealT, typename BufferT, typename IndexT>
inline void ReadNumpyArray(MatrixT<RealT>& m, BufferT* buf, IndexT shape0,
                           IndexT stride0, IndexT shape1, IndexT stride1) {
  auto m_it = m.reshaped().begin();
  for (IndexT i0 = 0; i0 < shape0; ++i0, buf += stride0) {
    BufferT* buf0 = buf;
    for (IndexT i1 = 0; i1 < shape1; ++i1, buf0 += stride1, ++m_it) {
      *m_it = *reinterpret_cast<RealT*>(buf0);
    }
  }
}
template <typename RealT, typename BufferT, typename IndexT>
inline void ReadNumpyArray(MatrixMapT<RealT>& m, BufferT* buf, IndexT shape0,
                           IndexT stride0, IndexT shape1, IndexT stride1) {
  new (&m) MatrixMapT<RealT>{reinterpret_cast<RealT*>(buf), shape1, shape0,
                             MatrixStrideT{stride0, stride1}};
}
template <typename RealT, typename BufferT, typename IndexT>
inline void ReadNumpyArray(VectorT<RealT>& m, BufferT* buf, IndexT shape0,
                           IndexT stride0) {
  auto m_it = m.reshaped().begin();
  for (IndexT i0 = 0; i0 < shape0; ++i0, buf += stride0, ++m_it) {
    *m_it = *reinterpret_cast<RealT*>(buf);
  }
}
template <typename RealT, typename BufferT, typename IndexT>
inline void ReadNumpyArray(VectorMapT<RealT>& m, BufferT* buf, IndexT shape0,
                           IndexT stride0) {
  new (&m) VectorMapT<RealT>{reinterpret_cast<RealT*>(buf), shape0,
                             VectorStrideT{stride0}};
}

/// Writes to buffer given contents of matrix/vector type
template <typename RealT, typename BufferT, typename IndexT>
inline void WriteNumpyArray(const MatrixT<RealT>& m, BufferT* buf,
                            IndexT shape0, IndexT stride0, IndexT shape1,
                            IndexT stride1) {
  auto m_it = m.reshaped().begin();
  for (IndexT i0 = 0; i0 < shape0; ++i0, buf += stride0) {
    BufferT* buf0 = buf;
    for (IndexT i1 = 0; i1 < shape1; ++i1, buf0 += stride1, ++m_it) {
      *reinterpret_cast<RealT*>(buf0) = *m_it;
    }
  }
}
template <typename RealT, typename BufferT, typename IndexT>
inline void WriteNumpyArray(const MatrixMapT<RealT>& m, BufferT* buf,
                            IndexT shape0, IndexT stride0, IndexT shape1,
                            IndexT stride1) {
  return;  // noop
}
template <typename RealT, typename BufferT, typename IndexT>
inline void WriteNumpyArray(const VectorT<RealT>& m, BufferT* buf,
                            IndexT shape0, IndexT stride0) {
  auto m_it = m.reshaped().begin();
  for (IndexT i0 = 0; i0 < shape0; ++i0, buf += stride0, ++m_it) {
    *reinterpret_cast<RealT*>(buf) = *m_it;
  }
}
template <typename RealT, typename BufferT, typename IndexT>
inline void WriteNumpyArray(const VectorMapT<RealT>& m, BufferT* buf,
                            IndexT shape0, IndexT stride0) {
  return;  // noop
}

template <typename RealT, NumpyStrideType t, typename BufferT, typename IndexT>
NumpyMatrixT<RealT, t> InitNumpyArray(BufferT* buf, IndexT shape0,
                                      IndexT stride0, IndexT shape1,
                                      IndexT stride1) {
  if constexpr (t == NumpyStrideType::CORE_ALIGNED) {
    return MatrixMapT<RealT>{reinterpret_cast<RealT*>(buf), shape1, shape0,
                             MatrixStrideT{stride0, stride1}};
  } else {
    MatrixT<RealT> result{shape1, shape0};
    if constexpr (t == NumpyStrideType::BROADCAST_ZERO) {
      ReadNumpyArray(result, buf, shape0, stride0, shape1, stride1);
    }
    return result;
  }
}
template <typename RealT, NumpyStrideType t, typename BufferT, typename IndexT>
NumpyVectorT<RealT, t> InitNumpyArray(BufferT* buf, IndexT shape0,
                                      IndexT stride0) {
  if constexpr (t == NumpyStrideType::CORE_ALIGNED) {
    return VectorMapT<RealT>{reinterpret_cast<RealT*>(buf), shape0,
                             VectorStrideT{stride0}};
  } else {
    VectorT<RealT> result{shape0};
    if constexpr (t == NumpyStrideType::BROADCAST_ZERO) {
      ReadNumpyArray(result, buf, shape0, stride0);
    }
    return result;
  }
}

}  // namespace moccasin
