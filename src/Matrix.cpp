//
// This file is part of the course materials for AMATH483/583 at the University
// of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0
// International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include "Matrix.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <thread>
using namespace std;

#ifdef __AVX__
#include "immintrin.h"
#endif    // __AVX__

void basicMultiply(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < B.num_cols(); ++j) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
}

void basicMultiplyMV(const Matrix& A, const Vector& x, Vector& y) {
  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      y(i) += A(i, j) * x(j);
    }
  }
}

void basicThreadedMultiply(const Matrix& A, const Matrix& B, Matrix& C) {
  size_t num_of_threads = 2;
  assert(A.num_rows() % num_of_threads == 0);
  std::vector<std::thread> threads(num_of_threads - 1);
  for (size_t thread = 1; thread < num_of_threads; ++thread) {
    threads[thread - 1] = std::thread([&, thread]() {
      for (size_t i = thread * A.num_rows() / num_of_threads; i < (thread + 1) * A.num_rows() / num_of_threads; ++i) {
        for (size_t j = 0; j < B.num_cols(); ++j) {
          for (size_t k = 0; k < A.num_cols(); ++k) {
            C(i, j) += A(i, k) * B(k, j);
          }
        }
      }
    });
  }

  for (size_t i = 0; i < A.num_rows() / num_of_threads; ++i) {
    for (size_t j = 0; j < B.num_cols(); ++j) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }

  for (size_t thread = 1; thread < num_of_threads; ++thread) {
    threads[thread - 1].join();
  }
}

void hoistedMultiply(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < B.num_cols(); ++j) {
      double t = C(i, j);
      for (size_t k = 0; k < A.num_cols(); ++k) {
        t += A(i, k) * B(k, j);
      }
      C(i, j) = t;
    }
  }
}

void tiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); i += 2) {
    for (size_t j = 0; j < B.num_cols(); j += 2) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
        C(i, j + 1) += A(i, k) * B(k, j + 1);
        C(i + 1, j) += A(i + 1, k) * B(k, j);
        C(i + 1, j + 1) += A(i + 1, k) * B(k, j + 1);
      }
    }
  }
}

void tiledMultiply2x4(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); i += 2) {
    for (size_t j = 0; j < B.num_cols(); j += 4) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
        C(i, j + 1) += A(i, k) * B(k, j + 1);
        C(i, j + 2) += A(i, k) * B(k, j + 2);
        C(i, j + 3) += A(i, k) * B(k, j + 3);
        C(i + 1, j) += A(i + 1, k) * B(k, j);
        C(i + 1, j + 1) += A(i + 1, k) * B(k, j + 1);
        C(i + 1, j + 2) += A(i + 1, k) * B(k, j + 2);
        C(i + 1, j + 3) += A(i + 1, k) * B(k, j + 3);
      }
    }
  }
}

void tiledMultiply4x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); i += 4) {
    for (size_t j = 0; j < B.num_cols(); j += 2) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
        C(i, j + 1) += A(i, k) * B(k, j + 1);
        C(i + 1, j) += A(i + 1, k) * B(k, j);
        C(i + 1, j + 1) += A(i + 1, k) * B(k, j + 1);
        C(i + 2, j) += A(i + 2, k) * B(k, j);
        C(i + 2, j + 1) += A(i + 2, k) * B(k, j + 1);
        C(i + 3, j) += A(i + 3, k) * B(k, j);
        C(i + 3, j + 1) += A(i + 3, k) * B(k, j + 1);
      }
    }
  }
}

void tiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); i += 4) {
    for (size_t j = 0; j < B.num_cols(); j += 4) {
      for (size_t k = 0; k < A.num_cols(); ++k) {
        C(i, j) += A(i, k) * B(k, j);
        C(i, j + 1) += A(i, k) * B(k, j + 1);
        C(i, j + 2) += A(i, k) * B(k, j + 2);
        C(i, j + 3) += A(i, k) * B(k, j + 3);
        C(i + 1, j) += A(i + 1, k) * B(k, j);
        C(i + 1, j + 1) += A(i + 1, k) * B(k, j + 1);
        C(i + 1, j + 2) += A(i + 1, k) * B(k, j + 2);
        C(i + 1, j + 3) += A(i + 1, k) * B(k, j + 3);
        C(i + 2, j) += A(i + 2, k) * B(k, j);
        C(i + 2, j + 1) += A(i + 2, k) * B(k, j + 1);
        C(i + 2, j + 2) += A(i + 2, k) * B(k, j + 2);
        C(i + 2, j + 3) += A(i + 2, k) * B(k, j + 3);
        C(i + 3, j) += A(i + 3, k) * B(k, j);
        C(i + 3, j + 1) += A(i + 3, k) * B(k, j + 1);
        C(i + 3, j + 2) += A(i + 3, k) * B(k, j + 2);
        C(i + 3, j + 3) += A(i + 3, k) * B(k, j + 3);
      }
    }
  }
}

void hoistedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (size_t i = 0; i < A.num_rows(); i += 2) {
    for (size_t j = 0; j < B.num_cols(); j += 2) {
      double t00 = C(i, j);
      double t01 = C(i, j + 1);
      double t10 = C(i + 1, j);
      double t11 = C(i + 1, j + 1);

      for (size_t k = 0; k < A.num_cols(); ++k) {
        t00 += A(i, k) * B(k, j);
        t01 += A(i, k) * B(k, j + 1);
        t10 += A(i + 1, k) * B(k, j);
        t11 += A(i + 1, k) * B(k, j + 1);
      }
      C(i, j)         = t00;
      C(i, j + 1)     = t01;
      C(i + 1, j)     = t10;
      C(i + 1, j + 1) = t11;
    }
  }
}

void blockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        for (size_t i = ii; i < ii + blocksize; i += 2) {
          for (size_t j = jj; j < jj + blocksize; j += 2) {
            for (size_t k = kk; k < kk + blocksize; ++k) {
              C(i, j) += A(i, k) * B(k, j);
              C(i, j + 1) += A(i, k) * B(k, j + 1);
              C(i + 1, j) += A(i + 1, k) * B(k, j);
              C(i + 1, j + 1) += A(i + 1, k) * B(k, j + 1);
            }
          }
        }
      }
    }
  }
}

void hoistedBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        for (size_t i = ii; i < ii + blocksize; i += 2) {
          for (size_t j = jj; j < jj + blocksize; j += 2) {

            double t00 = C(i, j);
            double t01 = C(i, j + 1);
            double t10 = C(i + 1, j);
            double t11 = C(i + 1, j + 1);

            for (size_t k = kk; k < kk + blocksize; ++k) {
              t00 += A(i, k) * B(k, j);
              t01 += A(i, k) * B(k, j + 1);
              t10 += A(i + 1, k) * B(k, j);
              t11 += A(i + 1, k) * B(k, j + 1);
            }

            C(i, j)         = t00;
            C(i, j + 1)     = t01;
            C(i + 1, j)     = t10;
            C(i + 1, j + 1) = t11;
          }
        }
      }
    }
  }
}

void copyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        Matrix BB(blocksize, blocksize);
        for (size_t j = jj, jb = 0; j < jj + blocksize; ++j, ++jb) {
          for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
            BB(jb, kb) = B(k, j);
          }
        }

        for (size_t i = ii; i < ii + blocksize; i += 2) {
          for (size_t j = jj, jb = 0; j < jj + blocksize; j += 2, jb += 2) {
            for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {

              C(i, j) += A(i, k) * BB(jb, kb);
              C(i, j + 1) += A(i, k) * BB(jb + 1, kb);
              C(i + 1, j) += A(i + 1, k) * BB(jb, kb);
              C(i + 1, j + 1) += A(i + 1, k) * BB(jb + 1, kb);
            }
          }
        }
      }
    }
  }
}

void hoistedCopyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        Matrix BB(blocksize, blocksize);
        for (size_t j = jj, jb = 0; j < jj + blocksize; ++j, ++jb) {
          for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
            BB(jb, kb) = B(k, j);
          }
        }

        for (size_t i = ii; i < ii + blocksize; i += 2) {
          for (size_t j = jj, jb = 0; j < jj + blocksize; j += 2, jb += 2) {

            double t00 = C(i, j);
            double t01 = C(i, j + 1);
            double t10 = C(i + 1, j);
            double t11 = C(i + 1, j + 1);

            for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
              t00 += A(i, k) * BB(jb, kb);
              t01 += A(i, k) * BB(jb + 1, kb);
              t10 += A(i + 1, k) * BB(jb, kb);
              t11 += A(i + 1, k) * BB(jb + 1, kb);
            }

            C(i, j)         = t00;
            C(i, j + 1)     = t01;
            C(i + 1, j)     = t10;
            C(i + 1, j + 1) = t11;
          }
        }
      }
    }
  }
}

#ifdef __AVX__
void hoistedCopyBlockedTiledMultiply2x2AVX(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        Matrix BB(blocksize, blocksize);
        for (size_t j = jj, jb = 0; j < jj + blocksize; ++j, ++jb) {
          for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
            BB(jb, kb) = B(k, j);
          }
        }

        for (size_t i = ii; i < ii + blocksize; i += 2) {
          for (size_t j = jj, jb = 0; j < jj + blocksize; j += 2, jb += 2) {

            __m256d tx = _mm256_setr_pd(C(i, j), C(i, j + 1), C(i + 1, j), C(i + 1, j + 1));

            for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
              __m256d a = _mm256_setr_pd(A(i, k), A(i, k), A(i + 1, k), A(i + 1, k));
              __m256d b = _mm256_setr_pd(BB(jb, kb), BB(jb + 1, kb), BB(jb, kb), BB(jb + 1, kb));
              __m256d c = _mm256_mul_pd(a, b);
              tx        = _mm256_add_pd(tx, c);
            }

            __m128d a = _mm256_extractf128_pd(tx, 0);
            __m128d b = _mm256_extractf128_pd(tx, 1);
            _mm_store_pd(&C(i, j), a);
            _mm_store_pd(&C(i + 1, j), b);
          }
        }
      }
    }
  }
}
#endif    // __AVX__

void hoistedCopyBlockedTiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        Matrix BB(blocksize, blocksize);
        for (size_t j = jj, jb = 0; j < jj + blocksize; ++j, ++jb) {
          for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
            BB(jb, kb) = B(k, j);
          }
        }

        for (size_t i = ii; i < ii + blocksize; i += 4) {
          for (size_t j = jj, jb = 0; j < jj + blocksize; j += 4, jb += 4) {

            double t00 = C(i, j);
            double t01 = C(i, j + 1);
            double t10 = C(i + 1, j);
            double t11 = C(i + 1, j + 1);
            double t20 = C(i + 2, j);
            double t21 = C(i + 2, j + 1);
            double t30 = C(i + 3, j);
            double t31 = C(i + 3, j + 1);
            double t02 = C(i, j + 2);
            double t03 = C(i, j + 3);
            double t12 = C(i + 1, j + 2);
            double t13 = C(i + 1, j + 3);
            double t22 = C(i + 2, j + 2);
            double t23 = C(i + 2, j + 3);
            double t32 = C(i + 3, j + 2);
            double t33 = C(i + 3, j + 3);

            for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
              t00 += A(i, k) * BB(jb, kb);
              t10 += A(i + 1, k) * BB(jb, kb);
              t20 += A(i + 2, k) * BB(jb, kb);
              t30 += A(i + 3, k) * BB(jb, kb);

              t01 += A(i, k) * BB(jb + 1, kb);
              t11 += A(i + 1, k) * BB(jb + 1, kb);
              t21 += A(i + 2, k) * BB(jb + 1, kb);
              t31 += A(i + 3, k) * BB(jb + 1, kb);

              t02 += A(i, k) * BB(jb + 2, kb);
              t12 += A(i + 1, k) * BB(jb + 2, kb);
              t22 += A(i + 2, k) * BB(jb + 2, kb);
              t32 += A(i + 3, k) * BB(jb + 2, kb);

              t03 += A(i, k) * BB(jb + 3, kb);
              t13 += A(i + 1, k) * BB(jb + 3, kb);
              t23 += A(i + 2, k) * BB(jb + 3, kb);
              t33 += A(i + 3, k) * BB(jb + 3, kb);
            }

            C(i, j)         = t00;
            C(i, j + 1)     = t01;
            C(i + 1, j)     = t10;
            C(i + 1, j + 1) = t11;
            C(i + 2, j)     = t20;
            C(i + 2, j + 1) = t21;
            C(i + 3, j)     = t30;
            C(i + 3, j + 1) = t31;
            C(i, j + 2)     = t02;
            C(i, j + 3)     = t03;
            C(i + 1, j + 2) = t12;
            C(i + 1, j + 3) = t13;
            C(i + 2, j + 2) = t22;
            C(i + 2, j + 3) = t23;
            C(i + 3, j + 2) = t32;
            C(i + 3, j + 3) = t33;
          }
        }
      }
    }
  }
}

#ifdef __AVX__
void hoistedCopyBlockedTiledMultiply4x4AVX(const Matrix& A, const Matrix& B, Matrix& C) {
  const size_t blocksize = std::min<size_t>(A.num_rows(), 64UL);

  for (size_t ii = 0; ii < A.num_rows(); ii += blocksize) {
    for (size_t jj = 0; jj < B.num_cols(); jj += blocksize) {
      for (size_t kk = 0; kk < A.num_cols(); kk += blocksize) {

        Matrix BB(blocksize, blocksize);
        for (size_t j = jj, jb = 0; j < jj + blocksize; ++j, ++jb) {
          for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {
            BB(jb, kb) = B(k, j);
          }
        }

        for (size_t i = ii; i < ii + blocksize; i += 4) {
          for (size_t j = jj, jb = 0; j < jj + blocksize; j += 4, jb += 4) {

            __m256d t0x = _mm256_load_pd(&C(i, j));
            __m256d t1x = _mm256_load_pd(&C(i + 1, j));
            __m256d t2x = _mm256_load_pd(&C(i + 2, j));
            __m256d t3x = _mm256_load_pd(&C(i + 3, j));

            for (size_t k = kk, kb = 0; k < kk + blocksize; ++k, ++kb) {

              __m256d bx = _mm256_setr_pd(BB(jb, kb), BB(jb + 1, kb), BB(jb + 2, kb), BB(jb + 3, kb));

              // Broadcast sholud be implemented with VBROADCASTSD
              __m256d a0 = _mm256_broadcast_sd(&A(i, k));
              a0         = _mm256_mul_pd(bx, a0);
              t0x        = _mm256_add_pd(t0x, a0);

              __m256d a1 = _mm256_broadcast_sd(&A(i + 1, k));
              a1         = _mm256_mul_pd(bx, a1);
              t1x        = _mm256_add_pd(t1x, a1);

              __m256d a2 = _mm256_broadcast_sd(&A(i + 2, k));
              a2         = _mm256_mul_pd(bx, a2);
              t2x        = _mm256_add_pd(t2x, a2);

              __m256d a3 = _mm256_broadcast_sd(&A(i + 3, k));
              a3         = _mm256_mul_pd(bx, a3);
              t3x        = _mm256_add_pd(t3x, a3);
            }

            _mm256_store_pd(&C(i, j), t0x);
            _mm256_store_pd(&C(i + 1, j), t1x);
            _mm256_store_pd(&C(i + 2, j), t2x);
            _mm256_store_pd(&C(i + 3, j), t3x);
          }
        }
      }
    }
  }
}
#endif    // __AVX__

Vector operator*(const Matrix& A, const Vector& x) {
  assert(A.num_cols() == x.num_rows());

  Vector y(A.num_rows(), 0.0);
  basicMultiplyMV(A, x, y);

  return y;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
  assert(A.num_cols() == B.num_rows());

  Matrix C(A.num_rows(), B.num_cols());
  zeroize(C);

  basicMultiply(A, B, C);

  return C;
}

Matrix operator+(const Matrix& A, const Matrix& B) {
  assert(A.num_cols() == B.num_cols() && A.num_rows() == B.num_rows());

  Matrix C(A.num_rows(), A.num_cols());

  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      C(i, j) = A(i, j) + B(i, j);
    }
  }
  return C;
}

Matrix operator-(const Matrix& A, const Matrix& B) {
  assert(A.num_cols() == B.num_cols() && A.num_rows() == B.num_rows());

  Matrix C(A.num_rows(), A.num_cols());

  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

Matrix transpose(const Matrix& A) {
  Matrix B(A.num_cols(), A.num_rows());
  for (size_t i = 0; i < B.num_rows(); ++i) {
    for (size_t j = 0; j < B.num_cols(); ++j) {
      B(i, j) = A(j, i);
    }
  }
  return B;
}

double frobeniusNorm(const Matrix& A) {
  double sum = 0.0;
  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      sum += A(i, j) * A(i, j);
    }
  }

  return sqrt(sum);
}

double oneNorm(const Matrix& A) {
  std::vector<double> v(A.num_cols(), 0.0);

  for (size_t j = 0; j < A.num_cols(); ++j) {
    for (size_t i = 0; i < A.num_rows(); ++i) {
      v[j] += std::abs(A(i, j));
    }
  }
  return *std::max_element(v.begin(), v.end());
}

double infinityNorm(const Matrix& A) {
  std::vector<double> v(A.num_rows(), 0.0);

  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      v[i] += std::abs(A(i, j));
    }
  }
  return *std::max_element(v.begin(), v.end());
}

void randomize(Matrix& A) {
  std::default_random_engine             generator;
  std::uniform_real_distribution<double> distribution(2.0, 32.0);
  static auto                            dice = std::bind(distribution, generator);

  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      A(i, j) = dice();
    }
  }
}

void zeroize(Matrix& C) {
  for (size_t i = 0; i < C.num_rows(); ++i) {
    for (size_t j = 0; j < C.num_cols(); ++j) {
      C(i, j) = 0.0;
    }
  }
}

void piscetize(Matrix& A, size_t xpoints, size_t ypoints) {
  assert(A.num_rows() == A.num_cols());
  assert(xpoints * ypoints == A.num_rows());

  for (size_t j = 0; j < xpoints; j++) {
    for (size_t k = 0; k < ypoints; k++) {
      size_t jrow = j * ypoints + k;

      if (j != 0) {
        size_t jcol   = (j - 1) * ypoints + k;
        A(jrow, jcol) = -1.0;
      }
      if (k != 0) {
        size_t jcol   = j * ypoints + (k - 1);
        A(jrow, jcol) = -1.0;
      }

      A(jrow, jrow) = 4.0;

      if (k != ypoints - 1) {
        size_t jcol   = j * ypoints + (k + 1);
        A(jrow, jcol) = -1.0;
      }
      if (j != xpoints - 1) {
        size_t jcol   = (j + 1) * ypoints + k;
        A(jrow, jcol) = -1.0;
      }
    }
  }
}

void writeMatrix(const Matrix& A, const string& filename) {
  ofstream outputFile(filename);
  streamMatrix(A, outputFile);
  outputFile.close();
}

void streamMatrix(const Matrix& A) { streamMatrix(A, cout); }

void streamMatrix(const Matrix& A, ostream& outputFile) {
  // Write header
  outputFile << "AMATH 583 MATRIX" << endl;
  outputFile << A.num_rows() << " " << A.num_cols() << endl;

  // Write data
  for (size_t i = 0; i < A.num_rows(); ++i) {
    for (size_t j = 0; j < A.num_cols(); ++j) {
      outputFile << A(i, j) << " ";
    }
    outputFile << endl;
  }

  // Write tailer
  outputFile << "THIS IS THE END" << endl;
}
