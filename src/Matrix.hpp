//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//
#if !defined(__MATRIX_HPP)
#define __MATRIX_HPP

#include "Vector.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

class Matrix {
public:
  explicit Matrix(size_t M, size_t N) : num_rows_(M), num_cols_(N), storage_(num_rows_ * num_cols_) {}
  explicit Matrix(size_t M, size_t N, double init) : num_rows_(M), num_cols_(N), storage_(num_rows_ * num_cols_, init) {}

  double& operator()(size_t i, size_t j) {
    assert(i < num_rows_ && j < num_cols_);

    return storage_[i * num_cols_ + j];
  }
  const double& operator()(size_t i, size_t j) const {
    assert(i < num_rows_ && j < num_cols_);
    return storage_[i * num_cols_ + j];
  }

  std::vector<double>::iterator       begin() { return storage_.begin(); }
  std::vector<double>::const_iterator begin() const { return storage_.begin(); }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return num_cols_; }

private:
  size_t              num_rows_, num_cols_;
  std::vector<double> storage_;
};

class MatrixView {
public:
  MatrixView(const Matrix& A, size_t beginrow, size_t endrow, size_t begincol, size_t endcol)
      : num_rows_(endrow - beginrow), num_cols_(endcol - begincol), array_iter(A.begin() + beginrow * A.num_cols() + begincol),
        lda(A.num_cols()) {}

  const double& operator()(size_t i, size_t j) const { return array_iter[i * lda + j]; }
  double&       operator()(size_t i, size_t j) { return const_cast<double&>(array_iter[i * lda + j]); }

  const size_t& num_rows() const { return const_cast<const size_t&>(num_rows_); }
  const size_t& num_cols() const { return const_cast<const size_t&>(num_cols_); }

protected:
  size_t                              num_rows_, num_cols_;
  std::vector<double>::const_iterator array_iter;
  size_t                              lda;
};

class TransposedMatrixView {
public:
  TransposedMatrixView(const Matrix& A, size_t beginrow, size_t endrow, size_t begincol, size_t endcol)
      : num_rows_(endcol - begincol), num_cols_(endrow - beginrow), array_iter(A.begin() + beginrow * A.num_cols() + begincol),
        lda(A.num_cols()) {}

  const double& operator()(size_t i, size_t j) const {
    assert(i < num_rows_ && j < num_cols_);
    return array_iter[j * lda + i];
  }

  double& operator()(size_t i, size_t j) {
    assert(i < num_rows_ && j < num_cols_);
    return const_cast<double&>(array_iter[j * lda + i]);
  }

  const size_t& num_rows() const { return const_cast<const size_t&>(num_rows_); }
  const size_t& num_cols() const { return const_cast<const size_t&>(num_cols_); }

protected:
  size_t                              num_rows_, num_cols_;
  std::vector<double>::const_iterator array_iter;
  size_t                              lda;
};

Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
double frobeniusNorm(const Matrix& A);
void   basicMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   basicThreadedMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedMultiply(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   blockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply2x4(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply4x2(const Matrix& A, const Matrix& B, Matrix& C);
void   tiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C);
void   copyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply2x2AVX(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C);
void   hoistedCopyBlockedTiledMultiply4x4AVX(const Matrix& A, const Matrix& B, Matrix& C);
double oneNorm(const Matrix& A);
double infinityNorm(const Matrix& A);
double frobeniusNorm(const Matrix& A);
void   zeroize(Matrix& C);
void   randomize(Matrix& A);
void   piscetize(Matrix& A, size_t xpoints, size_t ypoints);
void   writeMatrix(const Matrix& A, const std::string& filename);
void   streamMatrix(const Matrix& A);
void   streamMatrix(const Matrix& A, std::ostream& outputFile);

Vector operator*(const Matrix& A, const Vector& x);
void   basicMultiplyMV(const Matrix& A, const Vector& x, Vector& y);

#endif    // __MATRIX_HPP
