#pragma once
#include <cstddef>
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

namespace Isorropia {

template <typename T, size_t NumCols>
class mdspan {
private:
    T* data_;
    size_t num_rows_;

public:
    KOKKOS_INLINE_FUNCTION
    mdspan(T* ptr, size_t num_rows) : data_(ptr), num_rows_(num_rows) {}

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t row, size_t col) const {
        return data_[row * NumCols + col];
    }

    KOKKOS_INLINE_FUNCTION
    size_t num_rows() const { return num_rows_; }

    KOKKOS_INLINE_FUNCTION
    size_t num_cols() const { return NumCols; }
};

} // namespace Isorropia
