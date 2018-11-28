#include "Matrix.hpp"
#include <cstddef>
#include <future>

void matmat_iijjkk_64x64x64_BT_ij_H_k_2x2x1(const Matrix& A, const Matrix& B, Matrix& C) {

  size_t ii_begin = 0, ii_last = C.num_rows();
  size_t jj_begin = 0, jj_last = C.num_cols();
  size_t kk_begin = 0, kk_last = A.num_cols();

  size_t thread_count = 0;
  size_t num_threads = 8;
  std::vector<std::future<void>> futs (num_threads);

  for (size_t ii = ii_begin, iib = 0; ii < ii_last; ii += 64, iib += 64) {
    for (size_t jj = jj_begin, jjb = 0; jj < jj_last; jj += 64, jjb += 64) {

      futs[thread_count++] = std::async(std::launch::async, [&] () -> void {

      for (size_t kk = kk_begin, kkb = 0; kk < kk_last; kk += 64, kkb += 64) {

        size_t i_begin = ii, i_step = std::min<size_t>(64, ii_last - ii), i_end = i_begin + 2 * (i_step / 2),
               i_last  = i_begin + i_step;
        size_t j_begin = jj, j_step = std::min<size_t>(64, jj_last - jj), j_end = j_begin + 2 * (j_step / 2),
               j_last  = j_begin + j_step;
        size_t k_begin = kk, k_step = std::min<size_t>(64, kk_last - kk), k_end = k_begin + 1 * (k_step / 1),
               k_last = k_begin + k_step;

        Matrix B_T(j_step, k_step);
        for (size_t jb = 0, jt = j_begin; jb < j_step; ++jb, ++jt) {
          for (size_t kb = 0, kt = k_begin; kb < k_step; ++kb, ++kt) {
            B_T(jb, kb) = B(kt, jt);
          }
        }

        for (size_t i = i_begin, ib = 0; i < i_end; i += 2, ib += 2) {
          for (size_t j = j_begin, jb = 0; j < j_end; j += 2, jb += 2) {
            double t00 = C(i + 0, j + 0);
            double t01 = C(i + 0, j + 1);
            double t10 = C(i + 1, j + 0);
            double t11 = C(i + 1, j + 1);
            for (size_t k = k_begin, kb = 0; k < k_end; k += 1, kb += 1) {
              t00 += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              t01 += A(i + 0, k + 0) * B_T(jb + 1, kb + 0);
              t10 += A(i + 1, k + 0) * B_T(jb + 0, kb + 0);
              t11 += A(i + 1, k + 0) * B_T(jb + 1, kb + 0);
            }
            C(i + 0, j + 0) = t00;
            C(i + 0, j + 1) = t01;
            C(i + 1, j + 0) = t10;
            C(i + 1, j + 1) = t11;
            for (size_t k = k_end, kb = k_end - k_begin; k < k_last; ++k, ++kb) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 0, j + 1) += A(i + 0, k + 0) * B_T(jb + 1, kb + 0);
              C(i + 1, j + 0) += A(i + 1, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 1, j + 1) += A(i + 1, k + 0) * B_T(jb + 1, kb + 0);
            }
          }
          for (size_t j = j_end, jb = j_end - j_begin; j < j_last; ++j, ++jb) {
            for (size_t k = k_begin, kb = 0; k < k_end; k += 1, kb += 1) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 1, j + 0) += A(i + 1, k + 0) * B_T(jb + 0, kb + 0);
            }
            for (size_t k = k_end, kb = k_end - k_begin; k < k_last; ++k, ++kb) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 1, j + 0) += A(i + 1, k + 0) * B_T(jb + 0, kb + 0);
            }
          }
        }
        for (size_t i = i_end, ib = i_end - i_begin; i < i_last; ++i, ++ib) {
          for (size_t j = j_begin, jb = 0; j < j_end; j += 2, jb += 2) {
            for (size_t k = k_begin, kb = 0; k < k_end; k += 1, kb += 1) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 0, j + 1) += A(i + 0, k + 0) * B_T(jb + 1, kb + 0);
            }
            for (size_t k = k_end, kb = k_end - k_begin; k < k_last; ++k, ++kb) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
              C(i + 0, j + 1) += A(i + 0, k + 0) * B_T(jb + 1, kb + 0);
            }
          }
          for (size_t j = j_end, jb = j_end - j_begin; j < j_last; ++j, ++jb) {
            for (size_t k = k_begin, kb = 0; k < k_end; k += 1, kb += 1) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
            }
            for (size_t k = k_end, kb = k_end - k_begin; k < k_last; ++k, ++kb) {
              C(i + 0, j + 0) += A(i + 0, k + 0) * B_T(jb + 0, kb + 0);
            }
          }
        }
      }

	  });
      if (thread_count == num_threads) {
	for (size_t i = 0; i < thread_count; ++i) {
	  futs[i].get();
	}
	thread_count = 0;
      }
      for (size_t i = 0; i < thread_count; ++i) {
	futs[i].get();
      }
      thread_count = 0;
    }
    for (size_t i = 0; i < thread_count; ++i) {
      futs[i].get();
    }
    thread_count = 0;
  }
  for (size_t i = 0; i < thread_count; ++i) {
    futs[i].get();
  }
  thread_count = 0;
}
