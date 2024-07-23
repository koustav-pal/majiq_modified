/**
 * NumpyMaskToIndices.hpp
 *
 * Convert boolean mask to index vector
 */
#pragma once

#include <vector>

namespace moccasin {

template <typename ValueT, typename BufferT, typename IndexT>
inline void UpdateIndices(std::vector<ValueT>& indices, BufferT* buf,
                          IndexT shape0, IndexT stride0) {
  indices.resize(0);  // resize to zero given new buffer
  for (IndexT i0 = 0; i0 < shape0; ++i0, buf += stride0) {
    if (*reinterpret_cast<bool*>(buf)) {
      indices.push_back(static_cast<ValueT>(i0));
    }
  }
}

}  // namespace moccasin
