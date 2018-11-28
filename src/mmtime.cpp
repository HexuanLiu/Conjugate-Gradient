#include "Matrix.hpp"
#include "Timer.hpp"
#include <functional>
#include <iostream>
#include <string>
#include <thread>

#include "matmat.hpp"

int main(int argc, char* argv[]) {

  std::cout << "This machine has " << std::thread::hardware_concurrency() << " cores" << std::endl;

  size_t maxsize = (argc > 1) ? std::stod(argv[1]) : 32;

  std::cout << "N\tN*N\tTime\tFlops" << std::endl;

  for (size_t size = 64; size <= maxsize; size *= 2) {
    size_t numruns = 128ULL * 1048ULL * 1048ULL / (size * size * size) + 2;

    Matrix A(size, size), B(size, size), C(size, size);
    randomize(A);
    randomize(B);
    randomize(C);

    Timer T;
    T.start();
    for (size_t i = 0; i < numruns; ++i) {
      matmat(A, B, C, 4);
    }
    T.stop();
    std::cout << size << "\t" << size * size << "\t" << T.elapsed() << "\t"
              << 2.0 * 1.e3 * numruns * size * size * size / T.elapsed() << std::endl;
  }

  return 0;
}
