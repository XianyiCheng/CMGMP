#pragma once
#include <modus/common/eigen.hpp>

namespace modus
{

/**
 * This class applies dimensionality reduction using the concept of projective
 * equivalence from matroid theory. Concretely, we perform the following pre-
 * processing steps
 *  1. rescale hyperplane normals (rows) such that ‖n‖ = 1
 *  2. project A into rowspace coordinates, i.e. A' = AR, where R is an
 *     orthonormal basis for the row space
 *  3. remove non-unique rows, up to re-orientation
 *  4. remove zero rows
 *
 * We also guarantee that following statements hold post-preprocessing
 *  1. normals do not change orientation,
 *  2. rearrangement indices are in sorted order 
 *  3. there are no normals such that ϵ <= ‖n₁ - n₂‖ <= kϵ. Note, ϵ is
 *     dynamically scaled until this is true.
 *
 */
class MatroidPreprocessor {
 public:
  // Input matrix of matroid M[A].
  Eigen::MatrixXd A_;

  // Minimally projectively equivalent matrix A' such that M[A'] = M[A].
  Eigen::MatrixXd A_prime_;

  // The rhs projection matrix in A' = AR, before any reordering.
  Eigen::MatrixXd projection_;

  // Array of indices which rearrange A into A'.
  Eigen::VectorXi index_array_;

  // Hyperplanes with normals with one epsilon, i.e. ‖n₁ - n₂‖ <= ϵ, are
  // considered the same.
  double epsilon_;

  // Hyperplanes with distances ϵ <= ‖n₁ - n₂‖ <= kϵ, where k is the exclusion
  // factor, will trigger a rescaling of ϵ upwards. This is to ensure subsequent
  // computations are not ambiguous.
  double exclusion_factor_;

  // Verbose print outs.
  bool verbose_;

  void SetVerbose(bool verbose) { verbose_ = verbose; }

  // Return array of indices which rearranges A into A'.
  Eigen::VectorXi GetIndices() { return index_array_; }

  // Return the rhs projection matrix in A' = AR, before reordering.
  Eigen::MatrixXd GetProjectionMatrix() { return projection_; }

  // Return the minimum dimensional A' such that M[A'] = M[A].
  Eigen::MatrixXd GetMinimallyProjectivelyEquivalentMatrix() { return A_prime_; }

  void Preprocess(const Eigen::MatrixXd& A, double epsilon, 
                  double exclusion_factor=1e1);
};

}