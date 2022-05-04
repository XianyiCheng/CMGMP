#include <modus/modes/preprocess/matroid_preprocessor.hpp>
#include <modus/common/linear_algebra.hpp>
#include <iostream>


void modus::MatroidPreprocessor::Preprocess
  (const Eigen::MatrixXd& A, double epsilon, double exclusion_factor)
{
  // Store variables.
  A_ = A;
  A_prime_ = A;
  epsilon_ = epsilon;
  exclusion_factor_ = exclusion_factor;
  // Project into rowspace.
  size_t n = A_prime_.rows();
  orth(A_prime_.transpose(), epsilon_, projection_);
  A_prime_ = A_prime_ * projection_;
  // Normalize hyperplanes.
  for (size_t i = 0; i < n; i++) {
    MODUS_ASSERT(A_prime_.row(i).norm() > epsilon_);
    A_prime_.row(i).normalize();
  }
  // Create index array to keep track of reorderings.
  index_array_.resize(n);
  for (size_t i = 0; i < n; i++) {
    index_array_[i] = i;
  }
  // Remove duplicate and nonzero rows until kϵ < ‖n₁ - n₂‖.
  while (true) {
    // Remove zero rows.
    Eigen::VectorXi idx = NonzeroRows(A_prime_, epsilon_);
    A_prime_ = GetRows(A_prime_, idx);
    index_array_ = GetRows(index_array_, idx);
    // Get unique rows.
    idx = UniqueRowsReoriented(A_prime_, epsilon_);
    A_prime_ = GetRows(A_prime_, idx);
    index_array_ = GetRows(index_array_, idx);
    // Check that the normal exclusion zone is empty.
    bool empty = true;
    for (size_t i = 0; i < index_array_.size(); i++) {
      for (size_t j = i + 1; j < index_array_.size(); j++) {
        double d = (A_prime_.row(i) - A_prime_.row(j)).norm();
        MODUS_ASSERT(d > epsilon_);
        if (d <= exclusion_factor_ * epsilon_) {
          empty = false;
        }
      }
    }
    if (empty) { 
      break;
    }
    // MODUS_ASSERT_WARNING(empty, "Warning: increasing epsilon from %f to %f", epsilon_, 2*epsilon_);
    epsilon_ *= 2;
  }
  // Reorder rows in increasing index order.
  Eigen::VectorXi idx = LexicographicArgSortRows(index_array_);
  A_prime_ = GetRows(A_prime_, idx);
  index_array_ = GetRows(index_array_, idx);
  // Verbose summary.
  if (verbose_) {
    std::cout 
    << "Matroid Preprocessor" << std::endl
    << PrettyMatrix(A_, "A ") << std::endl
    << "eps: " << epsilon << " -> " << epsilon_ << std::endl
    << "idx: " << index_array_.transpose() << std::endl
    << PrettyMatrix(A_prime_, "A'") << std::endl;
  }
}