#include "Matrix.hpp"
#include "Timer.hpp"
#include <iostream>
#include <string>
#include <functional>

#include "matmat2.hpp"


int main(int argc, char* argv[]) {

  size_t M = (argc > 1) ? std::stod(argv[1]) : 32;
  size_t N = (argc > 2) ? std::stod(argv[2]) : 24;
  size_t K = (argc > 3) ? std::stod(argv[3]) : 48;

  Matrix A(M, K), B(K, N);
  randomize(A);
  randomize(B);

  Matrix C(M, N);
  matmat_jki_3x2x4(A, B, C);

  Matrix D(M, N);
  basicMultiply(A, B, D);
  
  double diff = frobeniusNorm(C-D);

  std::cout << frobeniusNorm(C) << " " << frobeniusNorm(D) << " diff = " << diff << std::endl;
#if 0
										    streamMatrix(A);
										    streamMatrix(B);
										    streamMatrix(C);
										    streamMatrix(D);
#endif
  return 0;
}
