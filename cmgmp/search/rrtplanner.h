#ifndef TREE_H
#define TREE_H
    #include "tree.h"
#endif

#ifndef PLANNINGWORLD_H
#define PLANNINGWORLD_H
    #include "../worlds/PlanningWorld.h"
#endif

#ifndef PLANNEROPTION_H
#define PLANNEROPTION_H
    #include "planner_option.h"
#endif

// #define METHOD_QPCC 0
// #define METHOD_LCP 1
// #define METHOD_ALL 2
// #define METHOD_NONE 3


// class RRTPlannerOptions{

// public:
//     int method = METHOD_QPCC;
//     int sampleSO3 = true;
//     double goal_biased_prob = 0.8;
//     int max_samples = 5;
//     Vector3d rotation_sample_axis;

//     RRTPlannerOptions(){};
//     RRTPlannerOptions(const RRTPlannerOptions& opts){
        
//         this->method = opts.method;
//         this->sampleSO3 = opts.sampleSO3;
//         this->goal_biased_prob = opts.goal_biased_prob;
//         this->max_samples = opts.max_samples;
//         this->rotation_sample_axis = opts.rotation_sample_axis;
//     }

// };

class RRTPlanner{

public:

    PlanningWorld* pw;
    Tree* T;
    RRTPlannerOptions* options;

    //bounds for sampling object configurations
    Vector3d X_lb;
    Vector3d X_ub;

    // object weight
    Vector6d f_w;

    // unit length for extend function
    double eps_trans;
    double eps_angle;

    // object surface
    std::vector<ContactPoint> object_surface;
    // std::vector<VectorXi> couple_idxes;

    RRTPlanner(PlanningWorld* pw_);

    void Initialize(Vector3d q_lb, Vector3d q_ub, Vector7d q_start, Vector6d f_w, std::vector<ContactPoint>& surface, 
    double angle_weight, double translation_weight, double eps_trans, double eps_angle);

    void SetInitialManpulatorConfig(const VectorXd& config);

    bool randomRelocateFingers(const VectorXd& mnp_config, Vector7d x,
                const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
                VectorXd& new_config);

    bool randomRelocateFingers_Quasidynamic(const VectorXd& mnp_config, Vector7d x, Vector6d v,
                const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
                VectorXd& new_config);

    void Extend_QPCC(int near_idx, Vector7d x_rand);
    void Extend_QPCC_Quasidynamic(int near_idx, Vector7d x_rand);


    /*
    void Extend_All(int near_idx, Vector7d x_rand);
    void Extend_LCP(int near_idx, Vector7d x_rand);
    
    void Extend_QPCC_BACKWARD(int near_idx, Vector7d x_rand);
    void Extend_None(int near_idx, Vector7d x_rand);
    */
    void Search(const RRTPlannerOptions& opts, Vector7d x_goal, double goal_thr, std::vector<int>* node_path, bool& ifsuccess, double& planningtime, bool ifquasidynamic=false);
    void VisualizePath(const std::vector<int>& node_path);

};

