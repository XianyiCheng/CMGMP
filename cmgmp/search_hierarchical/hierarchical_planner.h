#ifndef HTREE_H
#define HTREE_H
    #include "hierarchical_tree.h"
#endif

#ifndef PLANNINGWORLD_H
#define PLANNINGWORLD_H
    #include "../worlds/PlanningWorld.h"
#endif

#ifndef PLANNEROPTION_H
#define PLANNEROPTION_H
    #include "../search/planner_option.h"
#endif

class HPlanner{

public:

    PlanningWorld* pw;
    HTree* T;
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

    HPlanner(PlanningWorld* pw_);

    void Initialize(Vector3d q_lb, Vector3d q_ub, Vector7d q_start, Vector6d f_w, std::vector<ContactPoint>& surface, 
    double angle_weight, double translation_weight, double eps_trans, double eps_angle);

    bool randomRelocateFingers(const VectorXd& mnp_config, Vector7d x,
                const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
                VectorXd& new_config);

    bool randomRelocateFingers_Quasidynamic(const VectorXd& mnp_config, Vector7d x, Vector6d v,
                const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
                VectorXd& new_config);

    void Extend_ObjectPath(int near_idx, Vector7d x_rand);

    void Search(const RRTPlannerOptions& opts, Vector7d x_goal, double goal_thr, std::vector<int>* node_path, bool& ifsuccess, double& planningtime, bool ifquasidynamic=false);
    void VisualizePath(const std::vector<int>& node_path);

    bool GenerateFingerPlan(const std::vector<int> node_path, int& feasible_nodes, 
    std::vector<Vector7d>* object_configs, std::vector<VectorXd>* mnp_configs);

    bool ForwardCheck_Quasistatic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path);

    bool ForwardCheck_Quasidynamic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path);

    bool ForwardCheck_Dynamic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path);

};

