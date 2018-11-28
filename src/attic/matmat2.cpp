#include "Matrix.hpp"
#include <cstddef>

void matmat_jki_3x2x4(const Matrix& A, const Matrix& B, Matrix& C) {
  size_t i_begin = 0, i_last = C.num_rows(), i_end = 4 * (i_last / 4);
  size_t j_begin = 0, j_last = C.num_cols(), j_end = 3 * (j_last / 3);
  size_t k_begin = 0, k_last = A.num_cols(), k_end = 2 * (k_last / 2);
  for (size_t j = j_begin; j < j_end; j += 3) {
    for (size_t k = k_begin; k < k_end; k += 2) {
      for (size_t i = i_begin; i < i_end; i += 4) {
        // inner fringe = 0, middle_fringe = 0, outer_fringe = 0
        // ilimit = 4, jlimit = 3, klimit = 2
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 0) += A(i + 0, k + 1) * B(k + 1, j + 0);
        C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
        C(i + 0, j + 1) += A(i + 0, k + 1) * B(k + 1, j + 1);
        C(i + 0, j + 2) += A(i + 0, k + 0) * B(k + 0, j + 2);
        C(i + 0, j + 2) += A(i + 0, k + 1) * B(k + 1, j + 2);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 0) += A(i + 1, k + 1) * B(k + 1, j + 0);
        C(i + 1, j + 1) += A(i + 1, k + 0) * B(k + 0, j + 1);
        C(i + 1, j + 1) += A(i + 1, k + 1) * B(k + 1, j + 1);
        C(i + 1, j + 2) += A(i + 1, k + 0) * B(k + 0, j + 2);
        C(i + 1, j + 2) += A(i + 1, k + 1) * B(k + 1, j + 2);
        C(i + 2, j + 0) += A(i + 2, k + 0) * B(k + 0, j + 0);
        C(i + 2, j + 0) += A(i + 2, k + 1) * B(k + 1, j + 0);
        C(i + 2, j + 1) += A(i + 2, k + 0) * B(k + 0, j + 1);
        C(i + 2, j + 1) += A(i + 2, k + 1) * B(k + 1, j + 1);
        C(i + 2, j + 2) += A(i + 2, k + 0) * B(k + 0, j + 2);
        C(i + 2, j + 2) += A(i + 2, k + 1) * B(k + 1, j + 2);
        C(i + 3, j + 0) += A(i + 3, k + 0) * B(k + 0, j + 0);
        C(i + 3, j + 0) += A(i + 3, k + 1) * B(k + 1, j + 0);
        C(i + 3, j + 1) += A(i + 3, k + 0) * B(k + 0, j + 1);
        C(i + 3, j + 1) += A(i + 3, k + 1) * B(k + 1, j + 1);
        C(i + 3, j + 2) += A(i + 3, k + 0) * B(k + 0, j + 2);
        C(i + 3, j + 2) += A(i + 3, k + 1) * B(k + 1, j + 2);
      }
      for (size_t i = i_end; i < i_last; ++i) {
        // inner fringe = 1, middle_fringe = 0, outer_fringe = 0
        // ilimit = 1, jlimit = 3, klimit = 2
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 0) += A(i + 0, k + 1) * B(k + 1, j + 0);
        C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
        C(i + 0, j + 1) += A(i + 0, k + 1) * B(k + 1, j + 1);
        C(i + 0, j + 2) += A(i + 0, k + 0) * B(k + 0, j + 2);
        C(i + 0, j + 2) += A(i + 0, k + 1) * B(k + 1, j + 2);
      }
    }
    for (size_t k = k_end; k < k_last; ++k) {
      for (size_t i = i_begin; i < i_end; i += 4) {
        // inner fringe = 0, middle_fringe = 1, outer_fringe = 0
        // ilimit = 4, jlimit = 3, klimit = 1
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
        C(i + 0, j + 2) += A(i + 0, k + 0) * B(k + 0, j + 2);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 1) += A(i + 1, k + 0) * B(k + 0, j + 1);
        C(i + 1, j + 2) += A(i + 1, k + 0) * B(k + 0, j + 2);
        C(i + 2, j + 0) += A(i + 2, k + 0) * B(k + 0, j + 0);
        C(i + 2, j + 1) += A(i + 2, k + 0) * B(k + 0, j + 1);
        C(i + 2, j + 2) += A(i + 2, k + 0) * B(k + 0, j + 2);
        C(i + 3, j + 0) += A(i + 3, k + 0) * B(k + 0, j + 0);
        C(i + 3, j + 1) += A(i + 3, k + 0) * B(k + 0, j + 1);
        C(i + 3, j + 2) += A(i + 3, k + 0) * B(k + 0, j + 2);
      }
      for (size_t i = i_end; i < i_last; ++i) {
        // inner fringe = 1, middle_fringe = 1, outer_fringe = 0
        // ilimit = 1, jlimit = 3, klimit = 1
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
        C(i + 0, j + 2) += A(i + 0, k + 0) * B(k + 0, j + 2);
      }
    }
  }
  for (size_t j = j_end; j < j_last; ++j) {
    for (size_t k = k_begin; k < k_end; k += 2) {
      for (size_t i = i_begin; i < i_end; i += 4) {
        // inner fringe = 0, middle_fringe = 0, outer_fringe = 1
        // ilimit = 4, jlimit = 1, klimit = 2
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 0) += A(i + 0, k + 1) * B(k + 1, j + 0);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 0) += A(i + 1, k + 1) * B(k + 1, j + 0);
        C(i + 2, j + 0) += A(i + 2, k + 0) * B(k + 0, j + 0);
        C(i + 2, j + 0) += A(i + 2, k + 1) * B(k + 1, j + 0);
        C(i + 3, j + 0) += A(i + 3, k + 0) * B(k + 0, j + 0);
        C(i + 3, j + 0) += A(i + 3, k + 1) * B(k + 1, j + 0);
      }
      for (size_t i = i_end; i < i_last; ++i) {
        // inner fringe = 1, middle_fringe = 0, outer_fringe = 1
        // ilimit = 1, jlimit = 1, klimit = 2
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 0, j + 0) += A(i + 0, k + 1) * B(k + 1, j + 0);
      }
    }
    for (size_t k = k_end; k < k_last; ++k) {
      for (size_t i = i_begin; i < i_end; i += 4) {
        // inner fringe = 0, middle_fringe = 1, outer_fringe = 1
        // ilimit = 4, jlimit = 1, klimit = 1
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
        C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
        C(i + 2, j + 0) += A(i + 2, k + 0) * B(k + 0, j + 0);
        C(i + 3, j + 0) += A(i + 3, k + 0) * B(k + 0, j + 0);
      }
      for (size_t i = i_end; i < i_last; ++i) {
        // inner fringe = 1, middle_fringe = 1, outer_fringe = 1
        // ilimit = 1, jlimit = 1, klimit = 1
        C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
      }
    }
  }
}
