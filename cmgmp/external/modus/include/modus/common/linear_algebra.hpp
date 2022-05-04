#pragma once
#include <modus/common/eigen.hpp>
#include <vector>
#include <string>

namespace modus
{

/**
 * @brief Get the Sign Vector object
 * 
 * @tparam Derived0 
 * @param x 
 * @param eps 
 * @return std::string 
 */
template <class Derived0>
std::string GetSignVector(const Eigen::MatrixBase<Derived0>& x,
                          double eps);

/**
 * @brief Get the Sign Vector object
 * 
 * @tparam Derived0 
 * @tparam Derived1 
 * @param x 
 * @param N 
 * @param T 
 * @param eps 
 * @return std::string 
 */
template <class Derived0, class Derived1>
std::string GetSignVector(const Eigen::MatrixBase<Derived0>& x,
                          const Eigen::MatrixBase<Derived1>& N, 
                          const Eigen::MatrixBase<Derived1>& T, 
                          double eps);

/**
 * @brief Convenience function for pretty printing matrices.
 * 
 * @tparam Derived 
 * @param A 
 * @return std::string 
 */
template <class Derived>
std::string PrettyMatrix(const Eigen::MatrixBase<Derived>& A, 
                         const std::string& var_name="",
                         size_t precision=6);

/**
 * @brief Return the index vector of A sorted in a stable lexicographic order.
 *
 * @tparam Derived 
 * @param A 
 * @return Eigen::VectorXi 
 */
template <class Derived>
Eigen::VectorXi LexicographicArgSortRows(const Eigen::MatrixBase<Derived>& A);

/**
 * @brief Return the index vector of A sorted in a stable order.
 * 
 * @tparam Derived 
 * @param A 
 * @return Eigen::VectorXi 
 */
template <class Vector>
Eigen::VectorXi ArgSort(const Vector& A);

/**
 * @brief Return the unique rows. This function uses lexicographical_argsort.
 *
 * @tparam Derived 
 * @param A 
 * @param eps 
 * @return Eigen::VectorXi 
 */
template <class Derived>
Eigen::VectorXi UniqueRows(const Eigen::MatrixBase<Derived>& A, double eps);

/**
 * @brief Return the unique rows up to sign. This function uses lexicographical
 * argsort.
 *
 * @tparam Derived 
 * @param A 
 * @param eps 
 * @return Eigen::VectorXi 
 */
template <class Derived>
Eigen::VectorXi UniqueRowsReoriented(const Eigen::MatrixBase<Derived>& A, double eps);

/**
 * @brief Get the rows indexed by idx.
 * 
 * @tparam Derived 
 * @param A 
 * @param idx 
 * @param out 
 */
template <class Derived>
void GetRows(const Eigen::MatrixBase<Derived>& A, const Eigen::VectorXi& idx,
              typename Derived::PlainObject& out);

/**
 * @brief Returns the rows indexed by idx.
 * 
 * @tparam Derived 
 * @param A 
 * @param idx 
 * @return Derived::PlainObject 
 */
template <class Derived>
typename Derived::PlainObject GetRows(const Eigen::MatrixBase<Derived>& A, 
                                       const Eigen::VectorXi& idx);

/**
 * @brief Return the indices of the non-zero rows of A.
 * 
 * @tparam Derived 
 * @param A 
 * @param eps 
 * @return Eigen::VectorXi 
 */
template <class Derived>
Eigen::VectorXi NonzeroRows(const Eigen::MatrixBase<Derived>& A, double eps);

/**
 * @brief Return indices where A[i] == z, A ∈ Rⁿ.
 * 
 * @tparam Derived 
 * @param A 
 * @param z 
 * @param eps 
 * @return Eigen::MatrixXi 
 */
template <class Derived>
Eigen::VectorXi ArgWhereEqual(const Eigen::MatrixBase<Derived>& A,
                                typename Derived::Scalar z, double eps);

/**
 * @brief 
 * 
 * @param a 
 * @return Eigen::Matrix3d 
 */
Eigen::Matrix3d Adjoint(const Eigen::Vector3d& a);

/**
 * @brief Invert transform in-place.
 * 
 * @param R 
 * @param t 
 */
void InverseTransform(Eigen::Matrix3d& R, Eigen::Vector3d& t);

}

std::vector<int> lexographic_argsort(const Eigen::MatrixXd& H);

// TODO 
Eigen::VectorXi unique_rows(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                             double eps);
// TODO 
void get_cols(const Eigen::MatrixXd& A, const Eigen::VectorXi& idx, Eigen::MatrixXd& out);
Eigen::MatrixXd get_cols(const Eigen::MatrixXd& A, const Eigen::VectorXi& idx);

void image_basis(const Eigen::MatrixXd& A, Eigen::MatrixXd& image, double eps);
void kernel_basis(const Eigen::MatrixXd& A, Eigen::MatrixXd& null, double eps);
void image_and_kernel_bases(const Eigen::MatrixXd& A, Eigen::MatrixXd& image, 
                            Eigen::MatrixXd& kernel, double eps);

Eigen::MatrixXd image_basis(const Eigen::MatrixXd& A, double eps);
Eigen::MatrixXd kernel_basis(const Eigen::MatrixXd& A, double eps);

/**
 * @brief Compute an orthonormal basis using SVD.
 * 
 * @param A 
 * @param image 
 */
void orth(const Eigen::MatrixXd& A, double eps, Eigen::MatrixXd& image);
Eigen::MatrixXd orth(const Eigen::MatrixXd& A, double eps);

/**
 * @brief Compute the kernel using SVD.
 * 
 * @param A 
 * @param kernel 
 */
void null(const Eigen::MatrixXd& A, double eps, Eigen::MatrixXd& kernel);
Eigen::MatrixXd null(const Eigen::MatrixXd& A, double eps);

/**
 * @brief Partition and project H_s, int_pt onto nullspace of H_c, where H_c =
 * {h ∈ H : h ⋅ int_pt = 0} and H_s = H \ H_c. 
 *
 * @param H         n×d hyperplane matrix
 * @param d         n×1 hyperplane offsets
 * @param eps       tolerance
 * @param H_c       n_c×d contacting hyperplanes
 * @param H_s       n_s×d separating hyperplanes
 * @param d_s       n_s×1 separating hyperplanes offsets
 * @param int_pt    interior point
 * @param index_c   indices of contacting hyperplanes
 * @param index_s   indices of separating hyperplanes
 */
void partproj_nullspace(
    const Eigen::MatrixXd& H, const Eigen::VectorXd& d, double eps,
    Eigen::MatrixXd& H_c, Eigen::MatrixXd& H_s, Eigen::VectorXd& d_s, 
    Eigen::VectorXd& int_pt, Eigen::VectorXi& index_c, Eigen::VectorXi& index_s
    );

/**
 * @brief Project points onto their affinely independent subspace.
 * 
 * @param pts       d×n points
 * @param eps       tolerance
 * @param pts_aff   d'×n points in affine subspace
 */
void project_affine(const Eigen::MatrixXd& pts, double eps, 
                    Eigen::MatrixXd& pts_aff);


#include <modus/common/linear_algebra-impl.hpp>