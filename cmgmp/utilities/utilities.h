#pragma once
#ifndef _MATH_ULTILITIES_H_
#define _MATH_ULTILITIES_H_

#ifdef __cplusplus
    #include <cmath>
    #include <iostream>
#else
    #include <math.h>
    #include <stdio.h>
#endif

#include <memory> // smart pointers
#include <vector>

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense>


/////////////////////////////////////////////////////////////////////////
//                   types and static variables
/////////////////////////////////////////////////////////////////////////
using namespace Eigen;

typedef Eigen::Vector3d Vector3d;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix4d Matrix4d;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Quaterniond Quaterniond;
typedef Eigen::Matrix<double, 7, 1> Vector7d;

typedef Eigen::Vector3f Vector3f;
typedef Eigen::Matrix3f Matrix3f;
typedef Eigen::Matrix4f Matrix4f;
typedef Eigen::MatrixXf MatrixXf;
typedef Eigen::VectorXf VectorXf;
typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Quaternionf Quaternionf;
typedef Eigen::Quaterniond Quaterniond;

/////////////////////////////////////////////////////////////////////////
//                          iostream
/////////////////////////////////////////////////////////////////////////
void stream_array_in(std::ostream &st, double *array, int length);
void stream_array_in(std::ostream &st, float *array, int length);
void stream_array_in(std::ostream &st, int *array, int length);

/////////////////////////////////////////////////////////////////////////
//                          scalar
/////////////////////////////////////////////////////////////////////////
void truncate(double *ele, const double _min, const double _max);
// Generate a random number in (0, 1).
double randd();
// Generate a random integer in [0, N)
int randi(int N);

// set seed using time
int set_rand_seed();
/////////////////////////////////////////////////////////////////////////
//                          vector&array
/////////////////////////////////////////////////////////////////////////

void buf_insert(const double ele, const int size, double * buf);
void copyArray(const float *src, float *dest, int dim);
void copyArray(const double *src, double *dest, int dim);
void setArray(float *array, float value, int dim);
void truncate(float *array, float min, float max, int dim);
double vec_max(const double * vec, const int size);
double vec_min(const double * vec, const int size);
double vec_mean(const double * vec, const int size);
double vec_slope(const double * x, const double * y,const int size);
/**
 * Find in the vector the element with maximum abs value.
 *
 * @param[in]  vec   The vector.
 * @param[in]  size  Size of the vector.
 * @param      id    Id of the maximum element
 *
 * @return     The maximum element
 */
double vec_max_abs(const double * vec, const int size, int *id);
// numerical differentiation with low pass filter
// input x, calculate dx/dt
// s/(as+1),
double diff_LPF(const double xdold, const double xnew, const double xold, const double T,const double a);
void truncate6f(Vector6f *v, float min, float max);
void truncate6d(Vector6d *v, double min, double max);
void truncate6d(Vector6d *v, const Vector6d &min, const Vector6d &max);
void stream_array_in6f(std::ostream &st, const Vector6f &array);
void stream_array_in6d(std::ostream &st, const Vector6d &array);
void double2float(const double *array_in, float *array_out, int n);
void float2double(const float *array_in, double *array_out, int n);
/**
 * Finds an element in a std vector.
 *
 * @param[in]  vec   The vector
 * @param[in]  ele   The element
 *
 * @return     index of the found element; -1 if not found.
 */
int findInVector(std::vector<int> vec, int ele);
int findInEigenVector(const Eigen::VectorXi &vec, int ele);

/**
 * Finds multiple elements in a vector. Return the common elements among the
 * two vectors.
 *
 * @param[in]  vec   The vector
 * @param[in]  eles  The elements to be checked
 *
 * @return     The common elements (not their indices)
 */
std::vector<int> findInVector(std::vector<int> vec, std::vector<int> eles);


// lines, plucker coordinates
Vector6d getPluckerLine(const Vector3d &p, const Vector3d &n);
double reciprocalProduct(const Vector6d &line1, const Vector6d &line2);
// note: this distance is signed.
double distBTPluckerLines(const Vector6d &line1, const Vector6d &line2);
// rad
double angleBTPluckerLines(const Vector6d &line1, const Vector6d &line2);
double distPoint2PluckerLine(const Vector3d &p, const Vector6d &line);


/////////////////////////////////////////////////////////////////////////
//                          Matrices
/////////////////////////////////////////////////////////////////////////
MatrixXd pseudoInverse(const MatrixXd &a,
  double epsilon = std::numeric_limits<double>::epsilon());

// compute Reduced row echelon form of A. In place computation.
// return the rank of A. After the computation, the first rank rows of A are
// the row space of the input A; the rest of A are zeros.
// This is basically the Gaussian elimination.
//
// TOL: acceptable magnitude of a pivot
int rref(MatrixXd *A, double TOL = 1e-9);
/**
 * Gram-Schmidt procedure. Compute an unitary basis of the rows of A.
 *
 * @param      A     The input matrix. Will be modified.
 * @param[in]  TOL   Tol
 *
 * @return     rank of A
 */
int rowSpace(MatrixXd *A, double TOL = 1e-9);
int nullSpace(MatrixXd *A, MatrixXd *nullA, double TOL = 1e-9);
/////////////////////////////////////////////////////////////////////////
//                          Robotics
/////////////////////////////////////////////////////////////////////////
/*  Frames/spaces:
        W: world frame
        T: current tool frame
        So: set tool frame with offset
        Tf: transformed generalized space
    Quantities:
        SE3: 4x4 homogeneous coordinates
        se3: 6x1 twist coordinate of SE3
        spt: 6x1 special twist: 3x1 position, 3x1 exponential coordinate for rotation
        td: 6x1 time derivative of twist.
        v: 6x1 velocity, either spatial or body

        wrench: 6x1 wrench. Makes work with body velocity
*/

Matrix3d wedge(const Vector3d &v);
Matrix4d wedge6(const Vector6d &t);

//
// Transformations
//

Eigen::Matrix3d aa2mat(const double theta, const Eigen::Vector3d n);
Matrix3d quat2SO3(const Quaterniond &q);
Matrix3d quat2SO3(double qw, double qx, double qy, double qz);
Matrix3d so32SO3(const Vector3d &v);
Vector3d SO32so3(const Matrix3d &R);
void so32quat(const Vector3d &so3, double *q);
void SO32quat(const Matrix3d &SO3, double *q);
Matrix4d pose2SE3(const double *pose);
Matrix4d posemm2SE3(const double *pose);
Matrix4d se32SE3(const Vector6d &twist);
Matrix4d spt2SE3(const Vector6d &spt);
Matrix4d SE3Inv(const Matrix4d &SE3);
Vector6d SE32se3(const Matrix4d &SE3);
Vector6d SE32spt(const Matrix4d &SE3);
Matrix6d SE32Adj(const Matrix4d &SE3);
void SE32Pose(const Matrix4d &SE3, double *pose);
void SE32Posemm(const Matrix4d &SE3, double *pose);
Eigen::Matrix3f quat2m(const Eigen::Quaternionf &q);
Eigen::Matrix3d rotX(double angle_rad);
Eigen::Matrix3d rotY(double angle_rad);
Eigen::Matrix3d rotZ(double angle_rad);
/**
 * @brief      Gets the rotation matrix from z vector. The x and y axes are
 *             choosen arbitrarily.
 *
 * @param[in]  z     { Unit vector measured in world frame }
 *
 * @return     Rotation matrix from world to the rotated frame.
 */
Eigen::Matrix3d getRFromZ(const Eigen::Vector3d &z);

//
// Quaternions
//
Eigen::Quaternionf QuatMTimes(const Eigen::Quaternionf &q1,
  const Eigen::Quaternionf &q2);
Eigen::Quaterniond QuatMTimes(const Eigen::Quaterniond &q1,
    const Eigen::Quaterniond &q2);
float angBTquat(const Eigen::Quaternionf &q1, const Eigen::Quaternionf &q2);
double angBTquat(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2);

//
// SE(3)
//
class CartesianPose
{
public:
  CartesianPose();
  ~CartesianPose();
  /**
   * Construct the pose from a 1x7 vector.
   *
   * @param[in]  pose  Pose vector. [x y z qw qx qy qz]
   */
  CartesianPose(std::vector<double> pose);
  CartesianPose(double *pose);
  /**
   * Constructs the pose from an Eigen Matrix. T must be either:
   *  a 4x4 homogeneous matrix, or
   *  a 7x1 vector.
   *
   * @param[in]  T     The 4x4 homogeneous matrix or 7x1 vector.
   */
  CartesianPose(const Eigen::MatrixXd &T);
  CartesianPose(const Eigen::Isometry3d &iso);
  CartesianPose(const Eigen::Quaterniond &q, const Eigen::Vector3d &p);
  CartesianPose(const Eigen::Matrix3d &R, const Eigen::Vector3d &p);

  // other constructors
  CartesianPose(CartesianPose&& gp);
  CartesianPose& operator=(CartesianPose&& gp); // move assignment
  CartesianPose(const CartesianPose &gp); // copy
  CartesianPose& operator=(const CartesianPose& gp); // copy assignment

  static CartesianPose Identity();

  // types
  using Ptr = std::shared_ptr<CartesianPose>;
  using ConstPtr = std::shared_ptr<const CartesianPose>;

  // setters/getters
  void setRotationMatrix(const Eigen::Matrix3d &R);
  void setQuaternion(const Eigen::Quaterniond &q);
  void setQuaternion(const std::vector<double> &q);
  void setXYZ(const Eigen::Vector3d &p);
  void setXYZ(const std::vector<double> &p);
  void setX(double);
  void setY(double);
  void setZ(double);
  void scaleXYZ(double scale);
  Eigen::Matrix3d getRotationMatrix() const;
  Eigen::Quaterniond getQuaternion() const;
  Eigen::Vector3d getXYZ() const;
  double getX() const;
  double getY() const;
  double getZ() const;
  Eigen::Vector3d getXAxis() const;
  Eigen::Vector3d getYAxis() const;
  Eigen::Vector3d getZAxis() const;
  Eigen::Matrix4d getTransformMatrix() const;
  Eigen::Isometry3d getIsometry3d() const;
  std::vector<double> getVector() const;

  // operators
  CartesianPose operator*(const CartesianPose &pose) const;
  CartesianPose inv() const;

  // Transformations
  Eigen::Vector3d transformVec(const Eigen::Vector3d &v) const;
  Eigen::Vector3d transformPoint(const Eigen::Vector3d &p) const;
  Eigen::Quaterniond transformQuat(const Eigen::Quaterniond &q) const;

  // metric

  /**
   * Distance between this pose and another pose
   *
   * @param[in]  pose   Another pose
   * @param[in]  ratio  scaling of rotation to length
   *
   * @return     angle*ratio + length
   */
  double distBTPose(const CartesianPose & pose, double ratio = 1.0) const;

  // MISC
  void print() const;
  void printPose() const;
  std::string poseString() const;
// public:
//   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  double qw_;
  double qx_;
  double qy_;
  double qz_;
  Eigen::Vector3d *p_;
  Eigen::Matrix3d *R_;
};

//
// Vector operations
//
/**
 * find the (signed) angle from vector x to vector b.
 *
 * @param[in]  x            initial vector.
 * @param[in]  b            final vector.
 * @param[in]  z            if specified, it is the rotation axis. It is used
 *                          to specify positive direction of rotation.
 * @param[in]  nonnegative  when z is present, nonnegative will decide the
 *                          range of return angle between [-pi, pi] or [0, 2pi]
 *
 * @return     The angle. If @p z is not given, range is [0, pi]. If @p z is
 *             given, range is [-pi, pi] (@p nonnegative = false) or
 *             [0, 2pi] (@p nonnegative = true).
 */
double angBTVec(Eigen::Vector3d x, Eigen::Vector3d b,
    Eigen::Vector3d z = Eigen::Vector3d::Zero(), bool nonnegative = false);

//
// Others
//

// Return the 6x6 jacobian matrix mapping from spt time derivative
//  to body velocity.
// Jac * spt time derivative = body velocity
Matrix6d JacobianSpt2BodyV(const Matrix3d &R);


/**
 * Linear programming.
 *   min C'x
 *   s.t. Ax <= b
 *        Ae x == be
 *        xl <= x <= xu
 * Elements of xl and xu can be inf or NaN to indicate no constraint. Only
 * when isfinite() return true would the bound be considered.
 *
 * @param[in]  C             cost to be minimized
 * @param[in]  A             Inequality constraints, can be empty
 * @param[in]  b
 * @param[in]  Ae            Equality constraints, can be empty
 * @param[in]  be
 * @param[in]  xl            Lower bound of variable, can be empty.
 * @param[in]  xu            Upper bound of variable, can be empty.
 * @param      xs            Stores the solution
 * @param      optimal_cost  The optimal cost
 *
 * @return     True if the problem is feasible
 */
bool lp(const Eigen::VectorXd &C, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu, Eigen::VectorXd *xs,
    double *optimal_cost);

#endif // _MATH_ULTILITIES_H_
