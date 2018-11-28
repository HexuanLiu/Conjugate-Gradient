#include <cstddef>
#include "Matrix.hpp"

void matmat_ijk_2x2x1(const Matrix &A, const Matrix& B, Matrix&C) {
  size_t i_begin = 0, i_end = 2*(C.num_rows()/2), i_last = C.num_rows();
  size_t j_begin = 0, j_end = 2*(C.num_cols()/2), j_last = C.num_cols();
  size_t k_begin = 0, k_end = A.num_cols(), k_last = k_end;
  
  for (size_t i = i_begin; i < i_end; i += 2) {

    for (size_t j = j_begin; j < j_end; j += 2) {
      for (size_t k = k_begin; k < k_end; k += 1) {
	C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
	C(i + 0, j + 1) += A(i + 0, k + 0) * B(k + 0, j + 1);
	C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
	C(i + 1, j + 1) += A(i + 1, k + 0) * B(k + 0, j + 1);
      } // k
    } // j

    for (size_t j = j_end; j < j_last; ++j) { // j cleanup 
      for (size_t k = k_begin; k < k_end; ++k) {
	C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
	C(i + 1, j + 0) += A(i + 1, k + 0) * B(k + 0, j + 0);
      }
    }
  } // i

  for (size_t i = i_end; i < i_last; ++i) { 

    for (size_t j = j_begin; j < j_end; ++j) {
      for (size_t k = k_begin; k < k_end; ++k) {
	C(i + 0, j + 0) += A(i + 0, k + 0) * B(k + 0, j + 0);
      }
    }
    for (size_t j = j_end; j < j_last; ++j) {
      for (size_t k = k_begin; k < k_end; ++k) {
	C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
}
