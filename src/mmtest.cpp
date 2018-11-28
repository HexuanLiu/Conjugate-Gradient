#include "Matrix.hpp"
#include "Timer.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "matmat.hpp"

int main(int argc, char* argv[]) {

  size_t M = (argc > 1) ? std::stod(argv[1]) : 32;
  size_t N = (argc > 2) ? std::stod(argv[2]) : 24;
  size_t K = (argc > 3) ? std::stod(argv[3]) : 48;

  Matrix A(M, K), B(K, N);
  randomize(A);
  randomize(B);

  Matrix D(M, N);
  basicMultiply(A, B, D);

  Matrix C(M, N);

  matmat(A, B, C, 3);
  double diff = frobeniusNorm(C - D);

  if (std::abs(diff) > 1) {
    std::cout << frobeniusNorm(C) << " " << frobeniusNorm(D) << " diff = " << diff << std::endl;
    return -1;
  }

  return 0;
}
