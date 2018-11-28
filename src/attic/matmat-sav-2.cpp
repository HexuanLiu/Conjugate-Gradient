#include "Matrix.hpp"
#include <cstddef>

void matmat_jik_4x2x1(const Matrix& A, const Matrix& B, Matrix& C) {
  size_t i_begin = 0, i_last = C.num_rows(), i_end = 2 * (i_last / 2);
  size_t j_begin = 0, j_last = C.num_cols(), j_end = 4 * (j_last / 4);
  size_t k_begin = 0, k_last = A.num_cols(), k_end = 1 * (k_last / 1);
  for (size_t j = j_begin; j < j_end; j += 4) {
    for (size_t i = i_begin; i < i_end; i += 2) {
      for (size_t k = k_begin; k < k_end; k += 1) {
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
        C(i + 0, j + 2) += A(i + 0, k + 0) * B(k + 0, j + 2);
        C(i + 0, j + 3) += A(i + 0, k + 0) * B(k + 0, j + 3);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 1) += A(i + 1, k + 0) * B(k + 0, j + 1);
        C(i + 1, j + 2) += A(i + 1, k + 0) * B(k + 0, j + 2);
        C(i + 1, j + 3) += A(i + 1, k + 0) * B(k + 0, j + 3);
      }                                            // k
      for (size_t k = k_end; k < k_last; ++k) {    // k cleanup
      }
    }                                            // i
    for (size_t i = i_end; i < i_last; ++i) {    // i cleanup
      for (size_t k = k_begin; k < k_end; ++k) {
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
      }
    }
  }                                            // j
  for (size_t j = j_end; j < j_last; ++j) {    // j cleanup
    for (size_t i = i_begin; i < i_end; ++i) {
      for (size_t k = k_begin; k < k_end; ++k) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
    for (size_t i = i_end; i < i_last; ++i) {
      for (size_t k = k_begin; k < k_end; ++k) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
}
