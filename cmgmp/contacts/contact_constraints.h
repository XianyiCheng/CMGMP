#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

struct ContactPoint {
    Vector3d p;
    Vector3d n;
    double d = 0;

    ContactPoint(){}
    
    ContactPoint( const ContactPoint& pt){
        p = pt.p;
        n = pt.n;
        d = pt.d;
    }

    ContactPoint(Vector3d position, Vector3d normal){
        p = position;
        n = normal;
    }

    ContactPoint(Vector3d position, Vector3d normal, double distance){
        p = position;
        n = normal;
        d = distance;
    }

    bool operator==(const ContactPoint& pt) const
    {
        return (((p - pt.p).norm() < 1e-4) && ((n - pt.n).norm() < 1e-4));
    }
};


class FrictionConeLinearization {
// generate separation planes
// D x = 0 seperate the space of x

public:
    int number_of_sliding_planes;
    Matrix<double,Dynamic,6> D;
    Matrix<double,Dynamic,6> D_dual;
    Matrix<double,Dynamic,6> cones;
    Matrix<double,Dynamic,6> vectors;
    Matrix<int,Dynamic,Dynamic> cone_modes;
    Matrix<int,Dynamic,Dynamic> vector_modes;
    // std::vector<Vector6d> D;
    // std::vector<Vector6d> D_dual;
    // std::vector<Vector62d> cones;
    // std::vector<Vector6d> vectors;
    // std::vector<VectorXi> cone_modes;
    // std::vector<VectorXi> vector_modes;

    FrictionConeLinearization(int n);
    MatrixXd getVectorsbyMode(const VectorXi& mode);

};

class ContactConstraints {

public:
    //std::vector<ContactPoint> environment_contacts;
    //std::vector<ContactPoint> manipulator_contacts;
    FrictionConeLinearization* friction_cone;
    Eigen::Matrix<double, 3, 6> basis;
    
    ContactConstraints(int num_sliding_planes);
    //updateEnvironmentContacts(const std::vector<ContactPoint>& envs);
    //updateManipulatorContacts(const std::vector<ContactPoint>& mnps);
    void ModeConstraints(const std::vector<ContactPoint>&pts ,VectorXi mode, double mu, Vector6d f_external, MatrixXd* A, VectorXd* b, MatrixXd* G, VectorXd* h);
    void NormalVelocityConstraints(const std::vector<ContactPoint>&pts, MatrixXd* A_pr, VectorXd* b_pr);
    void TangentVelocityConstraints(const std::vector<ContactPoint>&pts, MatrixXd* T_pr, VectorXd* t_pr);

};

void mergeManipulatorandEnvironmentConstraints(const MatrixXd& A_mnp, const VectorXd& b_mnp, const MatrixXd& G_mnp, const VectorXd& h_mnp,
    const MatrixXd& A_env, const VectorXd& b_env, const MatrixXd& G_env, const VectorXd& h_env,
    MatrixXd* A_pr, VectorXd* b_pr, MatrixXd* G_pr, VectorXd* h_pr);

void mergeManipulatorandEnvironmentConstraints_relax(const MatrixXd& A_mnp, const VectorXd& b_mnp, const MatrixXd& G_mnp, const VectorXd& h_mnp,
    const MatrixXd& A_env, const VectorXd& b_env, const MatrixXd& G_env, const VectorXd& h_env,
    MatrixXd* A_pr, VectorXd* b_pr, MatrixXd* G_pr, VectorXd* h_pr);

void deleteZeroRows(const MatrixXd& A, const VectorXd& b, MatrixXd* A_pr, VectorXd* b_pr);
void mergeDependentRows(const MatrixXd& A, const VectorXd& b, MatrixXd* A_pr, VectorXd* b_pr);

void copy_points(const std::vector<ContactPoint>& pts, std::vector<ContactPoint>* pts_new);

void print_contacts(const std::vector<ContactPoint>& pts);
