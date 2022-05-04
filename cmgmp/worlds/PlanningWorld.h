// #include "world.h"
// #include "contact_constraints.h"
// #include "utilities.h"
#ifndef _WORLD_TEMPLATE
#define _WORLD_TEMPLATE
  #include "WorldTemplate.h"
#endif
#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "../contacts/contact_constraints.h"
#endif

// #ifndef MANIPULATOR_H
// #define MANIPULATOR_H
//     #include "manipulator.h"
// #endif
#ifndef MANIPULATORS_MANIPULATORTEMPLATE   
#define MANIPULATORS_MANIPULATORTEMPLATE   
    #include "../manipulators/ManipulatorTemplate.h"
#endif


// planning world: stores properties that are related to planning, provide functions for planning

class PlanningWorld{

public:
    WorldTemplate* world;
    ContactConstraints* cons;

    double mu_env = 0.8;
    double mu_mnp = 0.8;
    double charac_len = 1.0; // charateristic length = pi/object_length. Object scale = 1-5, charac_len = 1

    Matrix6d object_inertia = Matrix6d::Identity();

    PlanningWorld();

    PlanningWorld(WorldTemplate* w);

    Vector6d EnvironmentConstrainedVelocity(const Vector6d& v_desire, 
        const std::vector<ContactPoint>& envs, const VectorXi& env_mode);

    bool isForceBalance(const std::vector<ContactPoint>& mnps, 
        const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o);

    Vector6d QPVelocity(const Vector6d& v_desire, const std::vector<ContactPoint>& mnps, 
        const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o, double wt, double wa);

    Vector6d QPVelocity_Quasidynamic(const Vector6d& v_desire, const std::vector<ContactPoint>& mnps, 
        const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o, double wt, double wa,
        double h_time);

    bool ForwardIntegration(const Vector7d& x_start, const Vector7d& x_goal, 
        const VectorXd& mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode, const Vector6d& f_ext_w,
        double wt, double wa, std::vector<Vector7d>* path, bool ifquasidynamic = false);
    
    bool ForwardIntegration_ObjectPath(const Vector7d& x_start, const Vector7d& x_goal, 
        const std::vector<ContactPoint>& envs_, const VectorXi& env_mode,
        double wt, double wa, std::vector<Vector7d>* path);

    // void addManipulatorGraphics();

    // void updateManipulatorGraphicsTransformation(const VectorXd& config, const Vector7d& object_pose);

    // void VisualizeObjectPath(const std::vector<Vector7d>& path, const VectorXd& mnp_config);
    // void VisualizeObjectnPoints(const Vector7d& x, const std::vector<ContactPoint>& pts);
    // void ShowAllObjectPath(const std::vector<Vector7d>& path);
};

