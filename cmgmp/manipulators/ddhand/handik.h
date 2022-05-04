#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "two_link_kinematics.h"
#include <vector>
#include <cmath>

using namespace Eigen;

// double radians(double degrees){
//     return M_PI*degrees/180.0;
// }

// double degrees(double radians){
//     return 180.0*radians/M_PI;
// }

class HandIK
{
    private:
        double l0_, l1_, l2_;
        TwoLinkKinematics kinematics_;


    public:
        HandIK(double l0, double l1, double l2): 
            l0_(l0), 
            l1_(l1), 
            l2_(l2), 
            kinematics_(TwoLinkKinematics(l1,l2)){}
    
    std::tuple<double, double, double> get_lengths(){
        return {l0_, l1_, l2_};
    }
    
    void set_lengths(double l0, double l1, double l2)
    {
        l0_ = l0;
        l1_ = l1;
        l2_ = l2;
        return;
    }
    
    std::tuple<Vector2d, Vector2d> fk(Vector2d th1, Vector2d th2)
    {
        auto [l0, l1, l2] = get_lengths();

        Vector2d b1 = Vector2d(-l0/2, 0);
        Vector2d b2 = Vector2d(l0/2, 0);
        auto p1 = kinematics_.fk(th1);
        auto p2 = kinematics_.fk(th2);
        p2[0] = -p2[0];
        p1 = p1 + b1;
        p2 = p2 + b2;
        return {p1, p2};
    }

    std::tuple<Vector2d, Vector2d> ik(Vector2d p1, Vector2d p2)
    {
        auto [l0, l1, l2] = get_lengths();
        auto b1 = Vector2d(-l0/2, 0);
        auto b2 = Vector2d(l0/2, 0);
        p1 = p1-b1;
        p2 = p2-b2;
        p2[0] = -p2[0];
        
        auto [s11, s12, nsol1] = kinematics_.ik(p1);
        auto [s21, s22, nsol2] = kinematics_.ik(p2);
        
        return {s12, s22};
    }

    bool is_inverse_kinematics(Vector2d pp1, Vector2d pp2){

        auto [l0, l1, l2] = get_lengths();
        auto b1 = Vector2d(-l0/2, 0);
        auto b2 = Vector2d(l0/2, 0);
        Vector2d p1 = pp1-b1;
        Vector2d p2 = pp2-b2;
        p2[0] = -p2[0];

        if(pp1.norm() > 1e-8){
            auto [s11, s12, nsol1] = kinematics_.ik(p1);
            if (nsol1 == ZERO_SOLUTION)
                return false;
            if (!kinematics_.inJointLimits(s12))
                return false;
        }   

        if(pp2.norm() > 1e-8){
            auto [s21, s22, nsol2] = kinematics_.ik(p2);
            if (nsol2 == ZERO_SOLUTION)
                return false;
            if (!kinematics_.inJointLimits(s22))
                return false;
        }      
        
        return true;
        
    }
    
    Vector2d to_base_zero(Vector2d joints)
    {
        joints[0] += ((M_PI/180.0)*140);
        joints[1] -= ((M_PI/180.0)*6.05);
        joints[0] = -joints[0];
        return(joints);
    }
    Vector2d to_calib_zero(Vector2d joints)
    {
        joints[0] = -joints[0];
        joints[0] -= ((M_PI/180.0)*140);
        joints[1] += ((M_PI/180.0)*6.05);
        return(joints);
    }
};
