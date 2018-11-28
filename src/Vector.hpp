//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#if !defined(__VECTOR_HPP)
#define __VECTOR_HPP

#include <iostream>
#include <vector>

class Vector {
public:
  explicit Vector(size_t M) : iRows(M), storage_(iRows) {}
  explicit Vector(size_t M, double init) : iRows(M), storage_(iRows, init) {}

  double&       operator()(size_t i) { return storage_[i]; }
  const double& operator()(size_t i) const { return storage_[i]; }

  size_t num_rows() const { return storage_.size(); }

private:
  size_t              iRows;
  std::vector<double> storage_;
};

double oneNorm(const Vector& v);
double twoNorm(const Vector& v);
double infinityNorm(const Vector& v);
Vector operator*(const double& a, const Vector& x);
Vector operator+(const Vector& x, const Vector& y);
Vector operator-(const Vector& x, const Vector& y);
double dot(const Vector& x, const Vector& y);
void   zeroize(Vector& v);
void   randomize(Vector& v);

#endif    // __VECTOR_HPP
