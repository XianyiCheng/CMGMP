#include "cmgmp/utilities/utilities.h"
#include "cmgmp/utilities/sample.h"
#include "cmgmp/contacts/contact_kinematics.h"
#include "cmgmp/contacts/contact_constraints.h"
#include "cmgmp/utilities/eiquadprog.hpp"

const static double PI = 3.1415926;

void test_uniform_quaternion(){
    printf("test uniform quaternion: ");
    Quaterniond q = generate_unit_quaternion();
    std::cout << q.w() << ", " << q.x() << ", "<< q.y() << ", " << q.x() << std::endl;
}

void test_sample_position(){
    printf("test sample position: ");
    Vector3d ub(1,1,2);
    Vector3d lb(-1,0,-2);
    Vector3d p;
    p = sample_position(ub, lb);
    std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
}

void test_steer_position(){
    printf("test steer position: \n");
    Vector3d ub(1,1,2);
    Vector3d lb(-1,0,-2);
    Vector3d p0;
    Vector3d p1;
    p0 = sample_position(ub, lb);
    p1 = sample_position(ub, lb);
    Vector3d p = steer_position(p0, p1, 0.5);
    std::cout << "p0: "<< p0[0] << ", " << p0[1] << ", " << p0[2] << std::endl;
    std::cout << "p1: "<< p1[0] << ", " << p1[1] << ", " << p1[2] << std::endl;
    std::cout << "p: "<< p[0] << ", " << p[1] << ", " << p[2] << std::endl;
}

void test_steer_orientation(){
    printf("test steer orientation: \n");
    Quaterniond q0(1,0,0,0);
    std::cout << q0.normalized().toRotationMatrix() << std::endl;
    Quaterniond q1 = generate_unit_quaternion();
    Quaterniond q = steer_quaternion(q0,q1, PI/4);
    std::cout << "q0: " << q0.w() << ", " << q0.x() << ", "<< q0.y() << ", " << q0.x() << std::endl;
    std::cout << "q1: " << q1.w() << ", " << q1.x() << ", "<< q1.y() << ", " << q1.x() << std::endl;
    std::cout << "q: " << q.w() << ", " << q.x() << ", "<< q.y() << ", " << q.x() << std::endl;
    std::cout << "max angle: " << PI/4 << ", angle between q0, q: " << angBTquat(q0, q) << std::endl;
}

void test_contact_jacobian(){
    printf("test contact jacobian: \n");
    Vector3d p(1,2,3);
    Vector3d n(0,0,1);
    Matrix6d adgco = contact_jacobian(p,n);
    std::cout << adgco << std::endl;
}

void test_lp(){
    printf("test linear programming GLPK: \n");
    
    VectorXd C(2);
    C(0) = -0.6;
    C(1) = -0.5;
    MatrixXd A(2,2);
    A << 1,2,3,1;
    VectorXd b(2);
    b << 1,2;
    MatrixXd Ae;
    VectorXd be;
    VectorXd xl;
    VectorXd xu;
    VectorXd xs(2); 
    double optimal_cost;
    bool result;
    result = lp(C, A, b, Ae, be, xl, xu, &xs, &optimal_cost);
    std::cout << "Result: " << result << ", optimal cost: " << optimal_cost << ", solution: " << xs(0) << "," << xs(1) << std::endl;
    return;
}

void test_lp_2(){
    printf("test linear programming GLPK: \n");
    
    VectorXd C = VectorXd::Zero(6);
    
    MatrixXd A(6,6);
    A << -1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        0, 0,-1, 0, 0, 1,
        0,-1, 0, 0, 1, 0,
        -1, 0, 0,-1, 0, 0,
        0, 0, 0, 0, 0, 0;

    VectorXd b(6);
    b << 9.46499, 3.22555, -0.0992189, -0, -0, -0;

    MatrixXd G(10,6);
    G <<    0, 0, 1, 0, 0, 0,
   1, 1, 0.8,   0, 0, 0,
  -1, 1, 0.8,   0, 0, 0,
   1,-1, 0.8,   0, 0, 0,
  -1,-1, 0.8,   0, 0, 0,
   0, 0, 0, 0, 0, 1,
   0, 0, 0, 1, 1, 0.8,
   0, 0, 0,-1, 1, 0.8,
   0, 0, 0, 1,-1, 0.8,
   0, 0, 0,-1,-1, 0.8;

    VectorXd h = VectorXd::Zero(10);

    VectorXd xl;
    VectorXd xu;
    VectorXd xs(6); 
    double optimal_cost;
    bool result;
    result = lp(C, -G, -h, A, b, xl, xu, &xs, &optimal_cost);
    std::cout << "Result: " << result << ", optimal cost: " << optimal_cost << ", solution: " << xs.transpose() << std::endl;
    return;
}

void test_friction_cone(){
    FrictionConeLinearization fc(2);
    std::cout << "D: " << fc.D << std::endl;
    std::cout << "D_dual: " <<fc.D_dual << std::endl;
    std::cout << "cones: " <<fc.cones << std::endl;
    std::cout << "vecs: " <<fc.vectors << std::endl;
    std::cout << "cone modes: " <<fc.cone_modes << std::endl;
    std::cout << "vec modes: " <<fc.vector_modes << std::endl;
}

void test_mode_constraints(){

    ContactConstraints cons(2);

    std::vector<ContactPoint> envs;

    Vector3d nn(0,0,1);
    for (int i = 0; i< 4; i++){
        ContactPoint ee;
        ee.n = nn;
        envs.push_back(ee);
    }

    envs[0].p = Vector3d(0.5,0.5,-0.5);
    envs[1].p = Vector3d(0.5,-0.5,-0.5);
    envs[2].p = Vector3d(-0.5,0.5,-0.5);
    envs[3].p = Vector3d(-0.5,-0.5,-0.5);

    VectorXi mode;
    mode.resize(12);
    mode << 0,1,0,1,1,1,0,0,1,1,0,0;

    Vector6d f_ext;
    f_ext << 0,0,-10,0,0,0;

    MatrixXd A_env;
    MatrixXd G_env;
    VectorXd b_env;
    VectorXd h_env;

    cons.ModeConstraints(envs, mode, 0.6, f_ext, &A_env, &b_env, &G_env, &h_env);

    // std::cout << "A_env: " << A_env << std::endl;
    // std::cout << "b_env: " << b_env << std::endl; 
    // std::cout << "G_env: " << G_env << std::endl; 
    // std::cout << "h_env: " << h_env << std::endl;  

    // 
    std::vector<ContactPoint> mnps;

    {
        ContactPoint mm;
        mm.p = Vector3d(0,0,0.5);
        mm.n = Vector3d(0,0,-1);
        mnps.push_back(mm);
    }
    {
        ContactPoint mm;
        mm.p = Vector3d(0,-0.5,0);
        mm.n = Vector3d(0,1,0);
        mnps.push_back(mm);
    }

    for (const auto& pt: mnps){
        std::cout << "depth:" << pt.d << std::endl;
    }

    MatrixXd A_mnp;
    MatrixXd G_mnp;
    VectorXd b_mnp;
    VectorXd h_mnp;

    VectorXi mmode;
    mmode.resize(6);
    mmode << 0,0,0,0,0,0;

    cons.ModeConstraints(mnps, mmode, 0.8, f_ext, &A_mnp, &b_mnp, &G_mnp, &h_mnp);

    // std::cout << "A_mnp: " << A_mnp << std::endl;
    // std::cout << "b_mnp: " << b_mnp << std::endl; 
    // std::cout << "G_mnp: " << G_mnp << std::endl; 
    // std::cout << "h_mnp: " << h_mnp << std::endl;  

    MatrixXd A;
    MatrixXd G;
    VectorXd b;
    VectorXd h;

    mergeManipulatorandEnvironmentConstraints(A_mnp, b_mnp, G_mnp, h_mnp, A_env, b_env, G_env, h_env, &A, &b, &G, &h);

    // std::cout << "A: " << A << std::endl;
    // std::cout << "b: " << b << std::endl; 
    // std::cout << "G: " << G << std::endl; 
    // std::cout << "h: " << h << std::endl;  

    Vector6d v_goal;
    v_goal << 1,-0.5,0,0,0,0;

    int n_var = A.cols();

    MatrixXd P(n_var, n_var);
    VectorXd p(n_var);
    P.setIdentity();
    P.block(6,6,n_var-6, n_var-6) = 0.01*P.block(6,6,n_var-6, n_var-6);
    p.block(0,0,6,1) = -v_goal;

    VectorXd x(n_var);

    std::cout << "f: " << solve_quadprog(P, p, A.transpose(), -b,  G.transpose(), -h, x) << std::endl;
    std::cout << "x: \n" << x << std::endl;
    
}

void test_combination(){
    std::vector<int> h = {2,3,4,5};

    std::vector<std::vector<int>> p;

    combinations(h, 2, &p);

    for(auto& pi: p){
        for (auto i: pi){
            std::cout << i << ", "; 
        }
        std::cout << std::endl;
    }
}

int main(){
    // set_rand_seed();
    test_uniform_quaternion();
    test_sample_position();
    test_steer_position();
    test_steer_orientation();
    
    test_lp();
    test_lp_2();

    test_contact_jacobian();
    test_friction_cone();
    test_mode_constraints();
    test_combination();
}