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

#include "Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <random>

using namespace std;

double oneNorm(const Vector& v) {
  double sum = 0.0;
  for (size_t i = 0; i < v.num_rows(); ++i)
    sum += std::abs(v(i));
  return sum;
}

double twoNorm(const Vector& v) {
  double sum = 0.0;
  for (size_t i = 0; i < v.num_rows(); ++i)
    sum += v(i) * v(i);
  return std::sqrt(sum);
}

double infinityNorm(const Vector& v) {
  double d = 0.0;
  for (size_t i = 0; i < v.num_rows(); ++i)
    d = std::max(d, std::abs(v(i)));
  return d;
}

Vector operator*(const double& a, const Vector& x) {
  Vector y(x.num_rows());
  for (size_t i = 0; i < x.num_rows(); ++i) {
    y(i) = a * x(i);
  }
  return y;
}

Vector operator+(const Vector& x, const Vector& y) {
  assert(x.num_rows() == y.num_rows());

  Vector z(x.num_rows());
  for (size_t i = 0; i < x.num_rows(); ++i) {
    z(i) = x(i) + y(i);
  }
  return z;
}

Vector operator-(const Vector& x, const Vector& y) {
  assert(x.num_rows() == y.num_rows());

  Vector z(x.num_rows());
  for (size_t i = 0; i < x.num_rows(); ++i) {
    z(i) = x(i) - y(i);
  }
  return z;
}

double dot(const Vector& x, const Vector& y) {
  assert(x.num_rows() == y.num_rows());
  double sum = 0.0;
  for (size_t i = 0; i < x.num_rows(); ++i) {
    sum += x(i) * y(i);
  }
  return sum;
}

void zeroize(Vector& v) {
  for (size_t i = 0; i < v.num_rows(); ++i) {
    v(i) = 0.0;
  }
}

void randomize(Vector& v) {
  static std::default_random_engine             generator;
  static std::uniform_real_distribution<double> distribution(2.0, 32.0);
  static auto                                   dice = std::bind(distribution, generator);

  for (size_t i = 0; i < v.num_rows(); ++i) {
    v(i) = dice();
  }
}
