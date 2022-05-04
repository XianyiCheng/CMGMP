#include "cmgmp/contacts/contact_mode_enumeration.h"
#include "cmgmp/contacts/contact_constraints.h"
// #include <glog/logging.h>
// #include <modus/common/logging.hpp>

void test_8_contacts(){
    ContactConstraints cons(2);

    std::vector<ContactPoint> envs;

    {
    Vector3d nn(0,1,0);
    for (int i = 0; i< 4; i++){
        ContactPoint ee;
        ee.n = nn;
        envs.push_back(ee);
    }
    }
    envs[0].p = 2*Vector3d(0.5,-0.5,0.5);
    envs[1].p = 2*Vector3d(0.5,-0.5,-0.5);
    envs[2].p = 2*Vector3d(-0.5,-0.5,-0.5);
    envs[3].p = 2*Vector3d(-0.5,-0.5,0.5);


    {
    Vector3d nn(-1,0,0);
    for (int i = 0; i< 4; i++){
        ContactPoint ee;
        ee.n = nn;
        envs.push_back(ee);
    }
    }

    envs[4].p = 2*Vector3d(-0.5,0.5,0.5);
    envs[5].p = 2*Vector3d(-0.5,0.5,-0.5);
    envs[6].p = 2*Vector3d(-0.5,-0.5,-0.5);
    envs[7].p = 2*Vector3d(-0.5,-0.5,0.5);


    int n_pts = envs.size();
    MatrixXd A(n_pts, 6);
    VectorXd b(n_pts);
    MatrixXd T;
    VectorXd t;

    cons.NormalVelocityConstraints(envs, &A, &b);
    cons.TangentVelocityConstraints(envs, &T, &t);

    std::vector<VectorXi> cs_modes;
    std::vector<VectorXi> all_modes;

    cs_mode_enumeration(A, &cs_modes);
    std::cout << cs_modes.size() << std::endl;
    for (auto& m:cs_modes){
        std::cout << m.transpose() << std::endl;
    }

    // std::cout << A << std::endl;
    // Eigen::IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    // std::cout << A.format(CommaInitFmt) << std::endl;
    

    // all_mode_enumeration(A,T,&all_modes);

    // VectorXi cs_mode(8);
    // cs_mode << 1,1,0,0,0,0,1,1;
    // ss_mode_enumeration(A, T, cs_mode, &all_modes);

    // std::cout << all_modes.size() << std::endl;
    // for (auto& m:all_modes){
    //     std::cout << m.transpose() << std::endl;
    // }
}

void test_1_contacts(){
    ContactConstraints cons(2);

    std::vector<ContactPoint> envs;

    Vector3d nn(0,1,0);

    ContactPoint ee;
    ee.n = nn;
    ee.p = 2*Vector3d(0.5,-0.5,0.5);
    envs.push_back(ee);

    int n_pts = envs.size();
    MatrixXd A(n_pts, 6);
    VectorXd b(n_pts);
    MatrixXd T;
    VectorXd t;

    cons.NormalVelocityConstraints(envs, &A, &b);
    cons.TangentVelocityConstraints(envs, &T, &t);

    std::vector<VectorXi> cs_modes;
    std::vector<VectorXi> all_modes;

    all_mode_enumeration(A,T,&all_modes);
    std::cout << all_modes.size() << std::endl;
    for (auto& m:all_modes){
        std::cout << m.transpose() << std::endl;
    }

    // VectorXi cs_m(1);
    // cs_m[0]=0;
    // ss_mode_enumeration(A, T, cs_m, &all_modes);
    
    // for (auto& m:all_modes){
    //     std::cout << m.transpose() << std::endl;
    // }
}

int main(int argc, char **argv){
    test_8_contacts();
}