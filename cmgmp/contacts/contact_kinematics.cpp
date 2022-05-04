#include"contact_kinematics.h"
#include <math.h>

Matrix6d contact_jacobian(const Vector3d &position, const Vector3d &normal){
    Vector3d z(0,0,1); 
    Vector3d rotvec = z.cross(normal);
    Matrix3d R;

    if ((rotvec.norm() < 1e-6) && (signbit(z[2]) == signbit(normal[2]))){    
        R.setIdentity();
    } else if ((rotvec.norm() < 1e-6) && (signbit(z[2]) != signbit(normal[2]))){
        R << -1,0,0,0,1,0,0,0,-1;
    }
    else{
        rotvec = rotvec*(1/rotvec.norm());
        double rotangle = acos(z.dot(normal)/normal.norm());
        AngleAxisd aaxis(rotangle, rotvec);
        R = aaxis.toRotationMatrix();
    }
    // Matrix4d T = Matrix4d::Identity();
    // T.block<3, 3>(0, 0) = R;
    // T.block<3, 1>(0, 3) = position;
    // Matrix6d Adgoc = SE32Adj(T);
    Matrix4d T_inv = Matrix4d::Identity();
    T_inv.block<3, 3>(0, 0) = R.transpose();
    T_inv.block<3, 1>(0, 3) = -R.transpose()*position;
    Matrix6d Adgco = SE32Adj(T_inv);
    return Adgco;
}


Matrix4d pose2SE3(const Vector7d& x){
    //x: px, py, pz, qx, qy, qz, qw
    Matrix4d T;
    T.setIdentity();
    T.block(0,0,3,3) = quat2SO3(x(6), x(3), x(4), x(5));
    T.block(0,3,3,1) = x.block(0,0,3,1);

    return T;
}

Vector7d SE32pose(const Matrix4d& T){
    //x: px, py, pz, qx, qy, qz, qw
    Vector7d x;
    double q[4];
    x.block(0,0,3,1) = T.block(0,3,3,1);
    SO32quat(T.block(0,0,3,3), q);
    x(3) = q[1];
    x(4) = q[2];
    x(5) = q[3];
    x(6) = q[0];
    return x;
}

Vector6d compute_rbvel_body(const Vector7d& x, const Vector7d& x_goal){

    Matrix4d T = pose2SE3(x);
    Matrix4d T_goal = pose2SE3(x_goal);
    Matrix4d T_b = SE3Inv(T)*T_goal;

    Matrix3d R;
    R = T_b.block(0,0,3,3);
    AngleAxisd aa(R);
    
    Vector6d v_b;
    v_b.block(3,0,3,1) = aa.angle()*aa.axis();
    // std::cout << T_b << std::endl;
    v_b.block(0,0,3,1) = T_b.block(0,3,3,1);
    //Vector6d v_w = SE32Adj(T)*v_b;
    return v_b;
}

Vector6d compute_rbvel_spatial(const Vector7d& x, const Vector7d& x_goal){

    Matrix4d T = pose2SE3(x);
    Vector6d v_b = compute_rbvel_body(x, x_goal);
    Vector6d v_w = SE32Adj(T)*v_b;

    return v_w;
}
