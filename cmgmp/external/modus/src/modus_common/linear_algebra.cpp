#include <modus/common/linear_algebra.hpp>
#include <modus/modes/geometry/incidence_graph.hpp>
// #include <glog/logging.h>
#include <iostream>
#include <iomanip>


static int DEBUG=0;

using namespace modus;

bool lexographic_less_than(const std::pair<Eigen::VectorXd, int>& lhs,
               const std::pair<Eigen::VectorXd, int>& rhs) {
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

std::vector<int> lexographic_argsort(const Eigen::MatrixXd& H) {
  // Convert to std::vector.
  std::vector<std::pair<Eigen::VectorXd, int> > H_;
  for (int i = 0; i < H.rows(); i++) {
    H_.push_back(std::make_pair(H.row(i), i));
  }
  // Sort lexographically.
  std::sort(H_.begin(), H_.end(), &lexographic_less_than);
  if (DEBUG) {
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    Eigen::MatrixXd H_o(H_.size(), H.cols());
    for (int i = 0; i < H_.size(); i++) {
      H_o.row(i) = H_[i].first.transpose();
    }
    // DLOG(INFO) << "\nlexographical order" << "\n"
    //        << std::fixed << std::setprecision(6) << std::setfill(' ')
    //        << H_o.format(CleanFmt) << std::endl;
  }
  // Extract indices.
  std::vector<int> I;
  for (int i = 0; i < H_.size(); i++) {
    I.push_back(H_[i].second);
  }
  return I;
}

// Eigen::VectorXi unique_rows(const Eigen::MatrixXd& A, double eps) {
//     // std::cout << "A\n" << A << std::endl;
//     std::vector<int> idx_lex = lexographic_argsort(A);
//     int n_unique = 1;
//     for (int i = 1; i < idx_lex.size(); i++) {
//         int i0 = idx_lex[i-1];
//         int i1 = idx_lex[i];
//         if ((A.row(i0) - A.row(i1)).norm() > eps) {
//             n_unique += 1;
//         }
//     }
//     // std::cout << "num unique " << n_unique << std::endl;
//     Eigen::VectorXi idx_unique(n_unique);
//     idx_unique.setZero();
//     idx_unique[0] = idx_lex[0];
//     int k = 1;
//     for (int i = 1; i < idx_lex.size(); i++) {
//         int i0 = idx_lex[i-1];
//         int i1 = idx_lex[i];
//         if ((A.row(i0) - A.row(i1)).norm() > eps) {
//             // std::cout << k << std::endl;
//             // std::cout << idx_unique.transpose() << std::endl;
//             idx_unique[k++] = i1;
//         }
//     }
//     return idx_unique;
// }

Eigen::VectorXi unique_rows(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
              double eps) {
  Eigen::MatrixXd Ab(A.rows(), A.cols() + 1);
  Ab.leftCols(A.cols()) = A;
  Ab.rightCols(1) = b;
  return UniqueRows(Ab, eps);
}

// void get_rows(const Eigen::MatrixXd& A, 
//               const Eigen::VectorXi& idx, 
//               Eigen::MatrixXd& out)
// {
//     int n = idx.size();
//     int c = A.cols();
//     out.resize(n, c);
//     for (int i = 0; i < n; i++) {
//         out.row(i) = A.row(idx[i]);
//     }
// }

void get_cols(const Eigen::MatrixXd& A, 
        const Eigen::VectorXi& idx, 
        Eigen::MatrixXd& out)
{
  int n = idx.size();
  int r = A.rows();
  out.resize(r, n);
  for (int i = 0; i < n; i++) {
    out.col(i) = A.col(idx[i]);
  }
}

// Eigen::MatrixXd get_rows(const Eigen::MatrixXd& A, const Eigen::VectorXi& idx) {
//     Eigen::MatrixXd out;
//     get_rows(A, idx, out);
//     return out;
// }

Eigen::MatrixXd get_cols(const Eigen::MatrixXd& A, const Eigen::VectorXi& idx) {
  Eigen::MatrixXd out;
  get_cols(A, idx, out);
  return out;
}

// void get_rows(const Eigen::VectorXd& A, 
//               const Eigen::VectorXi& idx, 
//               Eigen::VectorXd& out)
// {
//     int n = idx.size();
//     out.resize(n);
//     for (int i = 0; i < n; i++) {
//         out[i] = A[idx[i]];
//     }
// }

// Eigen::VectorXd get_rows(const Eigen::VectorXd& A, const Eigen::VectorXi& idx) {
//     Eigen::VectorXd out;
//     get_rows(A, idx, out);
//     return out;
// }

void image_basis(const Eigen::MatrixXd& A, Eigen::MatrixXd& image, double eps) {
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
  cod.setThreshold(eps);
  cod.compute(A);
  image = cod.matrixQ();
  // Check determinant = 1 to see if the column space basis is in SO(n). If
  // not, permute any two rows to flip the determinant.
  if (std::abs(image.determinant() - 1) >= eps) {
    image.col(0).swap(image.col(1));
    if (DEBUG) {
      assert(std::abs(image.determinant() - 1) < eps);
    }
  }
  image = image.leftCols(cod.rank());
}

void kernel_basis(const Eigen::MatrixXd& A, Eigen::MatrixXd& kernel, double eps) {
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
  cod.setThreshold(eps);
  cod.compute(A);
  Eigen::MatrixXd Z = cod.matrixZ();
  Eigen::MatrixXd P_inv = cod.colsPermutation().inverse();
  kernel = (Z * P_inv).bottomRows(cod.dimensionOfKernel()).transpose();
  // if (DEBUG) {
  if (false) {
    Eigen::MatrixXd Q = cod.matrixQ();
    Eigen::MatrixXd T = cod.matrixT().topLeftCorner(cod.rank(), cod.rank()).triangularView<Eigen::Upper>();
    std::cout << "A\n" << A << std::endl;
    std::cout << "Q\n" << Q << std::endl;
    std::cout << "T\n" << T << std::endl;
    std::cout << "ZP\n" << Z*P_inv << std::endl;
    std::cout << "TZP\n" << T*Z*P_inv << std::endl;
    std::cout << "QTZP\n" << Q*T*Z*P_inv << std::endl;
    std::cout << "kernel\n" << kernel << std::endl;
    std::cout << "A*(ZP)^T\n" << A*(Z*P_inv).transpose() << std::endl;
  }
  MODUS_ASSERT(kernel.rows() == A.cols(),
    "Error: computed kernel with incorrect dimensions\n%s\n  eps: %f", 
    PrettyMatrix(A, "    A:", 17).c_str(), eps);
}

void image_and_kernel_bases(const Eigen::MatrixXd& A, Eigen::MatrixXd& image,
              Eigen::MatrixXd& kernel, double eps) {
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
  cod.setThreshold(eps);
  cod.compute(A);
  image = cod.matrixQ();
  image = image.leftCols(cod.rank());
  Eigen::MatrixXd Z = cod.matrixZ();
  Eigen::MatrixXd P_inv = cod.colsPermutation().inverse();
  kernel = (Z * P_inv).bottomRows(cod.dimensionOfKernel()).transpose();
}

Eigen::MatrixXd image_basis(const Eigen::MatrixXd& A, double eps) {
  Eigen::MatrixXd image;
  image_basis(A, image, eps);
  return image;
}

Eigen::MatrixXd kernel_basis(const Eigen::MatrixXd& A, double eps) {
  Eigen::MatrixXd kernel;
  kernel_basis(A, kernel, eps);
  return kernel;
}

void orth(const Eigen::MatrixXd& A, double eps, Eigen::MatrixXd& image) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd;
  svd.setThreshold(eps);
  svd.compute(A, Eigen::ComputeFullU);
  // std::cout << "svd U\n" << svd.matrixU() << std::endl;
  // std::cout << "svd U det " << svd.matrixU().determinant() << std::endl;
  image = svd.matrixU().leftCols(svd.rank());
}

Eigen::MatrixXd orth(const Eigen::MatrixXd& A, double eps) {
  Eigen::MatrixXd image;
  orth(A, eps, image);
  return image;
}

void partproj_nullspace(const Eigen::MatrixXd& H, const Eigen::VectorXd& d, double eps,
  Eigen::MatrixXd& H_c, Eigen::MatrixXd& H_s, Eigen::VectorXd& d_s,
  Eigen::VectorXd& int_pt, Eigen::VectorXi& index_c, Eigen::VectorXi& index_s) {
  // Find hyperplanes which are always in contact.
  Eigen::VectorXd residual = H * int_pt;
  arg_where_zero(residual, eps, index_c);
  arg_where_nonzero(residual, eps, index_s);
  // Parition hyperplanes.
  GetRows(H, index_c, H_c);
  GetRows(H, index_s, H_s);
  GetRows(d, index_s, d_s);
  // Get nullspace. 
  Eigen::MatrixXd kernel;
  kernel_basis(H_c, kernel, eps);
  // Project hyperplanes and interior point into contacting nullspace.
  H_s = H_s * kernel; // Project hyperplanes.
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
  qr.setThreshold(eps);
  qr.compute(kernel);
  Eigen::VectorXd tmp = qr.solve(int_pt); // Project interior point.
  int_pt = tmp;

  // DLOG(INFO) << "\n" << "partproj nullspace"
  //            << "\n\t" << "projected residual: " << residual.transpose()
  //            << "\n\t" << "     index contact: " << index_c.transpose()
  //            << "\n\t" << "    index separate: " << index_s.transpose();
}

void project_affine(const Eigen::MatrixXd& pts, double eps, 
          Eigen::MatrixXd& pts_aff) {

  // Check if 0 ∈ Aff(pts)
  int d = pts.rows();
  int n = pts.cols();
  Eigen::MatrixXd A(d+1,n);
  A.topRows(d) = pts;
  A.bottomRows(1) = Eigen::VectorXd::Ones(n).transpose();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(d+1);
  b[d] = 1;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
  qr.setThreshold(1e-8);
  qr.compute(A);
  Eigen::VectorXd x = qr.solve(b);

  // If 0 ∈ Aff(pts) continue, else shift affine subspace to the origin.
  Eigen::VectorXd x0 = (A * x).topRows(d);
  double s = A.row(d) * x;
  if (x0.norm() < eps && std::abs(s - 1) < eps) {
    A = pts;
  } else {
    A = pts;
    Eigen::VectorXd p0 = A.col(0);
    for (int i = 0; i < pts.cols(); i++) {
      A.col(i) -= p0;
    }
  }

  // Project points onto coordinates of the image (column) space.
  Eigen::MatrixXd image = orth(A, eps);
  // image_basis(A, image, eps);
  pts_aff = image.transpose() * A;

  // Eigen::MatrixXd image2 = orth(A, eps);

  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "\t\t[", "]");
  // DLOG(INFO) << "\n" << "project affine" << std::fixed << std::setprecision(6) << std::setfill(' ')
  //            << "\n\t" << "pts\n" << pts.format(CleanFmt)
  //            << "\n\t" << "pts shift\n" << A.format(CleanFmt)
  //            << "\n\t" << "colspace\n" << image.format(CleanFmt)
  //         //    << "\n\t" << "det(col)\n\t\t" << image.determinant()
  //            << "\n\t" << "col.col\n" << (image.transpose() * image).format(CleanFmt)
  //            << "\n\t" << "pts affine\n" << pts_aff.format(CleanFmt);
  if (DEBUG) {
    // Assert 0 ∈ affine hull of projected points.
    int d_ = pts_aff.rows();
    Eigen::MatrixXd A_(d_+1, n);
    A_.topRows(d_) = pts_aff;
    A_.bottomRows(1) = Eigen::VectorXd::Ones(n).transpose();
    Eigen::VectorXd b_(d_+1);
    b_.setZero();
    b_[d_] = 1;
    qr.compute(A_);
    Eigen::VectorXd x_ = qr.solve(b_);
    // std::cout << A_ << std::endl;
    // std::cout << x_ << std::endl;
    // std::cout << A_*x_ << std::endl;
    assert((A_*x_).topRows(d_).norm() < eps);
    assert(std::abs(A_.row(d_)*x_-1) < eps);
    // Assert unprojected points match input.
    // std::cout << "\tmax error " << (A - image * pts_aff).cwiseAbs().maxCoeff() << std::endl;
    assert((A - image * pts_aff).norm() < eps);
  }

  // Eigen::MatrixXd orth(6,4);
  // orth << 
  // -0.,       0.,      -0.,       0.,     
  //  0.5688,   0.40368, -0.11623, -0.,     
  // -0.5688,  -0.40368,  0.11623, -0.,     
  //  0.44169, -0.75974, -0.47718, -0.,     
  // -0.3973,   0.31124, -0.8633,  -0.,     
  // -0.,       0.,       0.,      -1.;
  // std::cout << "orth\n" << orth.format(CleanFmt) << std::endl;
  // pts_aff = orth.transpose() * A;
  // std::cout << "pts_aff\n" << pts_aff.format(CleanFmt) << std::endl;

  // std::cout << "pts\n" << A.format(CleanFmt) << std::endl;
  // std::cout << "image\n" << image.format(CleanFmt) << std::endl;
  // std::cout << "pts aff\n" << pts_aff.format(CleanFmt) << std::endl;
  // DLOG(INFO) << "\n\t" << "     rank: " << qr.rank()
  //            << "\n\t" << "residuals: " << (A*x).transpose()
  //            << "\n\t" << " colspace: " << image.rows() << "x" << image.cols();
  // DLOG(INFO) << "\n\tpoints affine:\n" << pts_aff;
}

Eigen::Matrix3d modus::Adjoint(const Eigen::Vector3d& a)
{
  Eigen::Matrix3d A;
  A << 0, -a[2], a[1],
     a[2], 0, -a[0],
     -a[1], a[0], 0;
  return A;
}

void modus::InverseTransform(Eigen::Matrix3d& R, Eigen::Vector3d& t)
{
  R.transposeInPlace();
  t = -R*t;
}

void arg_where_zero(const Eigen::VectorXd& x, double eps, Eigen::VectorXi& idx) {
  int cnt = 0;
  for (int i = 0; i < x.size(); i++) {
    cnt += abs(x[i]) <= eps;
  }
  idx.resize(cnt);
  int k = 0;
  for (int i = 0; i < x.size(); i++) {
    if (abs(x[i]) <= eps) {
      idx[k++] = i;
    }
  }
}

void arg_where_nonzero(const Eigen::VectorXd& x, double eps, Eigen::VectorXi& idx) {
  int cnt = 0;
  for (int i = 0; i < x.size(); i++) {
    cnt += abs(x[i]) > eps;
  }
  idx.resize(cnt);
  int k = 0;
  for (int i = 0; i < x.size(); i++) {
    if (abs(x[i]) > eps) {
      idx[k++] = i;
    }
  }
}