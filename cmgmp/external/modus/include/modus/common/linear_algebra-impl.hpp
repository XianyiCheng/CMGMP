#pragma once
#include <modus/common/linear_algebra.hpp>
#include <modus/common/eigen.hpp>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <modus/common/assert.hpp>


template <class Derived0>
std::string modus::GetSignVector(const Eigen::MatrixBase<Derived0>& x,
                 double eps)
{
  MODUS_ASSERT(eps > 0);
  std::string sign_vector;
  for (size_t i = 0; i < x.size(); i++) {
    if (x[i] < -eps) {
      sign_vector.push_back('-');
      MODUS_ASSERT(x[i] < -10*eps, "Value within 1e1 of tolerance: %f", x[i]);
    } else if (x[i] > eps) {
      sign_vector.push_back('+');
      MODUS_ASSERT(x[i] > 10*eps, "Value within 1e1 of tolerance: %f", x[i]);
    } else {
      sign_vector.push_back('0');
    }
  }
  return sign_vector;
}

template <class Derived0, class Derived1>
std::string modus::GetSignVector(const Eigen::MatrixBase<Derived0>& x,
                 const Eigen::MatrixBase<Derived1>& N, 
                 const Eigen::MatrixBase<Derived1>& T, 
                 double eps)
{
  MODUS_ASSERT(T.rows() % N.rows() == 0);
  std::string cs_mode = modus::GetSignVector(N * x, eps);
  std::string ss_mode = modus::GetSignVector(T * x, eps);
  std::string mode = cs_mode;
  size_t n_t = T.rows() / N.rows();
  for (size_t i = 0; i < cs_mode.size(); i++) {
    if (cs_mode[i] == '0') {
      for (size_t j = 0; j < n_t; j++) {
        mode.push_back(ss_mode[i * n_t + j]);
      }
    }
  }
  return mode;
}

template <class Derived>
std::string modus::PrettyMatrix(const Eigen::MatrixBase<Derived>& A, 
                const std::string& var_name,
                size_t precision) {
  var_name.size();
  std::string prefix = var_name;
  // std::string prefix(' ', var_name.size());
  if (var_name.size() > 0) {
    prefix.push_back(' ');
  }
  prefix.push_back('[');
  Eigen::IOFormat CleanFmt(precision, 0, ", ", "\n", prefix, "]");
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << std::setfill(' ') 
     << std::showpos << A.format(CleanFmt);
  return ss.str();
}

template <class Scalar>
bool lexicographic_less_than(const std::pair<std::vector<Scalar>, int>& lhs,
               const std::pair<std::vector<Scalar>, int>& rhs) {
  assert(lhs.first.size() == rhs.first.size());
  for (int i = 0; i < lhs.first.size(); i++) {
    double l = lhs.first[i];
    double r = rhs.first[i];
    // TODO Add arguement for threshold
    if (std::abs(l-r) < 1e-8) {
      continue;
    } else if (l < r) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

template <class Derived>
Eigen::VectorXi modus::LexicographicArgSortRows(const Eigen::MatrixBase<Derived>& A) {
  // Convert into vector of row vectors.
  std::vector<std::pair<std::vector<typename Derived::Scalar>, int> > Avv;
  for (int i = 0; i < A.rows(); i++) {
    std::vector<typename Derived::Scalar> row;
    for (int j = 0; j < A.cols(); j++) {
      row.push_back(A(i,j));
    }
    Avv.push_back(std::make_pair(row, i));
  }
  // Sort lexicographically.
  std::sort(Avv.begin(), Avv.end(), &lexicographic_less_than<typename Derived::Scalar>);
  // Extract indices.
  Eigen::VectorXi idx(Avv.size());
  for (int i = 0; i < Avv.size(); i++) {
    idx[i] = Avv[i].second;
  }
  return idx;
}

template <class Scalar>
bool argsort_less_than(const std::pair<Scalar, int>& lhs,
                       const std::pair<Scalar, int>& rhs) 
{
  return lhs.first < rhs.first;
}

template <class Vector>
Eigen::VectorXi modus::ArgSort(const Vector& A) {
  // Create an identical vector with index meta-data.
  using Scalar = typename Vector::value_type;
  using ScalarIndex = std::pair<Scalar, int>;
  std::vector<ScalarIndex> A_;
  for (size_t i = 0; i < A.size(); i++) {
    A_.push_back(std::make_pair(A[i], i));
  }
  // Sort.
  std::sort(A_.begin(), A_.end(), &argsort_less_than<Scalar>);
  // Extract indices.
  Eigen::VectorXi idx(A_.size());
  for (size_t i = 0; i < A_.size(); i++) {
    idx[i] = A_[i].second;
  }
  return idx;
}

template <class Derived>
Eigen::VectorXi modus::UniqueRows(const Eigen::MatrixBase<Derived>& A, double eps) {
  // Count the number of unique elements.
  Eigen::VectorXi idx_lex = LexicographicArgSortRows(A);
  int n_unique = 1;
  for (int i = 1; i < idx_lex.size(); i++) {
    int i0 = idx_lex[i-1];
    int i1 = idx_lex[i];
    if ((A.row(i0) - A.row(i1)).norm() > eps) {
      n_unique += 1;
    }
  }
  // Get the indices of the first unique elements.
  Eigen::VectorXi idx_unique(n_unique);
  idx_unique.setZero();
  idx_unique[0] = idx_lex[0];
  int k = 1;
  for (int i = 1; i < idx_lex.size(); i++) {
    int i0 = idx_lex[i-1];
    int i1 = idx_lex[i];
    if ((A.row(i0) - A.row(i1)).norm() > eps) {
      idx_unique[k++] = i1;
    }
  }
  return idx_unique;
}

template <class Derived>
Eigen::VectorXi modus::UniqueRowsReoriented(const Eigen::MatrixBase<Derived>& A, double eps) {
  // Reorient the matrix such that all first elements are positive.
  Derived A_ = A;
  for (int i = 0; i < A.rows(); i++) {
    if (A_(i,0) < 0) {
      A_.row(i) *= -1;
    }
  }
  return UniqueRows(A_, eps);
}

template <class Derived>
void modus::GetRows(const Eigen::MatrixBase<Derived>& A, const Eigen::VectorXi& idx,
        typename Derived::PlainObject& out) {
  // 
  int n = idx.size();
  int c = A.cols();
  out.resize(n, c);
  for (int i = 0; i < n; i++) {
    out.row(i) = A.row(idx[i]);
  }
}

template <class Derived>
typename Derived::PlainObject modus::GetRows(const Eigen::MatrixBase<Derived>& A, 
                     const Eigen::VectorXi& idx) {
  typename Derived::PlainObject out;
  GetRows(A, idx, out);
  return out;
}

template <class Derived>
Eigen::VectorXi modus::NonzeroRows(const Eigen::MatrixBase<Derived>& A, double eps) {
  std::vector<int> idx_std;
  for(int i = 0; i < A.rows(); i++) {
    if (A.row(i).norm() > eps) {
      idx_std.push_back(i);
    }
  }
  Eigen::VectorXi idx(idx_std.size());
  for (int i = 0; i < idx_std.size(); i++) {
    idx[i] = idx_std[i];
  }
  return idx;
}

template <class Derived>
Eigen::VectorXi modus::ArgWhereEqual(const Eigen::MatrixBase<Derived>& A,
                typename Derived::Scalar z, double eps) {
  // 
  std::vector<int> idx_std;
  for (int i = 0; i < A.rows(); i++) {
    if (std::abs(A(i) - z) < eps) {
      idx_std.push_back(i);
    }
  }
  Eigen::VectorXi idx(idx_std.size());
  for (int i = 0; i < idx_std.size(); i++) {
    idx[i] = idx_std[i];
  }
  return idx;
}