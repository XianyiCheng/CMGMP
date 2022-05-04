#include "hierarchical_planner.h"

#ifndef CONTACTMODE_H
#define CONTACTMODE_H
    #include "../contacts/contact_mode_enumeration.h"
#endif

#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "../contacts/contact_constraints.h"
#endif

#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "../contacts/contact_kinematics.h"
#endif

#ifndef SAMPLE_H
#define SAMPLE_H
    #include "../utilities/sample.h"
#endif

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#include <ctime>

double dist_vel(const Vector6d& v, const Vector6d& v0, double wt, double wa){
    double d = wt*(v-v0).block(0,0,3,1).norm() + wa*(v-v0).block(3,0,3,1).norm();
    return d;
}

Vector6d weight_w2o(const Vector7d& x, const Vector6d& f_ext_w){
    Matrix4d T = pose2SE3(x);
    Matrix6d Adg = SE32Adj(T);

    Matrix4d T_;
    T_.setIdentity();
    T_.block(0,0,3,3) = T.block(0,0,3,3);
    Vector6d f_ext_o = SE32Adj(T_).transpose()*f_ext_w;
    return f_ext_o;
}

void copy_pts(const std::vector<ContactPoint>& pts, std::vector<ContactPoint>* pts_new){
    for(auto& pt:pts){
        pts_new->push_back(pt);
    }
}

Vector7d steer_config(Vector7d x_near, Vector7d x_rand, double epsilon_translation, double epsilon_angle){
    
    // double epsilon_translation = 1.5;
    // double epsilon_angle = 3.14*100/180;

    Vector3d p_rand = x_rand.head(3);
    Vector3d p_near = x_near.head(3);
    Quaterniond q_near(x_near[6], x_near[3], x_near[4], x_near[5]); 
    Quaterniond q_rand(x_rand[6], x_rand[3], x_rand[4], x_rand[5]); 
    p_rand = steer_position(p_near, p_rand, epsilon_translation);
    q_rand = steer_quaternion(q_near, q_rand, epsilon_angle);
    Vector7d x_steer;
    x_steer << p_rand[0], p_rand[1], p_rand[2], double(q_rand.x()),  double(q_rand.y()),  double(q_rand.z()),  double(q_rand.w());
    return x_steer;
}

HPlanner::HPlanner(PlanningWorld* pw_){
    this->pw = pw_;
}

void HPlanner::Initialize(Vector3d q_lb, Vector3d q_ub, Vector7d q_start, Vector6d f_w, std::vector<ContactPoint>& surface, 
    double angle_weight, double translation_weight, double eps_trans, double eps_angle)
{
    // angle_weight and translation_weight: heuristics of how important the dimensions are. 

    // this->T = new Tree(angle_weight, translation_weight);
    this->T = new HTree(angle_weight, this->pw->charac_len*translation_weight);
    this->X_ub = q_ub;
    this->X_lb = q_lb;
    this->f_w = f_w;
    this->eps_trans = eps_trans;
    this->eps_angle = eps_angle;

    this->object_surface = surface;

    HNode start_node(q_start);
    this->pw->world->getObjectContacts(&(start_node.envs), start_node.config);
    T->initial_node(&start_node); 

}

bool HPlanner::randomRelocateFingers(const VectorXd& mnp_config, Vector7d x,
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
    VectorXd& new_config)
{

    this->pw->world->updateObjectPose(x);

    VectorXd new_mnp_config(0);
    int n_pts = this->pw->world->getRobot()->getNumberofFingertips();

    VectorXi relocate_mode(env_mode.size());
    relocate_mode.setZero();
    
    int samples = 20;

    // for (int n_on = 1; n_on <= n_pts; n_on++){
        //n_on: number of relocating fingers
        for(int k_sample = 0; k_sample < samples; k_sample++){

            int n_on = randi(n_pts+1);
            
            std::vector<ContactPoint> remain_mnps;
            bool isresample = this->pw->world->getRobot()->resampleFingers(n_on, mnp_config, x, this->object_surface, new_mnp_config, &remain_mnps);
            
            if (!isresample){
                continue;
            }

            std::vector<ContactPoint> remain_grasp_mnps;
            this->pw->world->getRobot()->Fingertips2PointContacts(remain_mnps, &remain_grasp_mnps);
            bool isremainbalance = this->pw->isForceBalance(remain_grasp_mnps, envs, relocate_mode, f_o);
            
            std::vector<ContactPoint> grasp_mnps;
            std::vector<ContactPoint> mnp_fingertips;
            this->pw->world->getRobot()->getFingertipsOnObject(new_mnp_config, x, &mnp_fingertips);
            this->pw->world->getRobot()->Fingertips2PointContacts(mnp_fingertips, &grasp_mnps);
            
            bool isbalance = this->pw->isForceBalance(grasp_mnps, envs, env_mode, f_o);
            bool ifCollide = this->pw->world->isRobotCollide(new_mnp_config);

            if (isremainbalance && isbalance && !ifCollide){
                new_config = new_mnp_config;
                return true;
            } 
        }
    // }

    return false;
}

bool HPlanner::randomRelocateFingers_Quasidynamic(const VectorXd& mnp_config, Vector7d x, Vector6d v,
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
    VectorXd& new_config)
{

    this->pw->world->updateObjectPose(x);
    VectorXd new_mnp_config(0);
    int n_pts = this->pw->world->getRobot()->getNumberofFingertips();

    VectorXi relocate_mode(env_mode.size());
    relocate_mode.setZero();
    
    int samples = 10;

    // for (int n_on = 1; n_on <= n_pts; n_on++){
        //n_on: number of relocating fingers
        for(int k_sample = 0; k_sample < samples; k_sample++){

            int n_on = randi(n_pts+1);
            
            std::vector<ContactPoint> remain_mnps;
            bool isresample = this->pw->world->getRobot()->resampleFingers(n_on, mnp_config, x, this->object_surface, new_mnp_config, &remain_mnps);
            
            if (!isresample){
                continue;
            }

            std::vector<ContactPoint> remain_grasp_mnps;
            this->pw->world->getRobot()->Fingertips2PointContacts(remain_mnps, &remain_grasp_mnps);
            bool isremainbalance = this->pw->isForceBalance(remain_grasp_mnps, envs, relocate_mode, f_o);

            Vector6d v_dyn_remain = this->pw->QPVelocity_Quasidynamic(v, remain_grasp_mnps, envs,
                    env_mode, f_o, T->translation_weight, T->angle_weight,1);

            bool isremaindynamic = v.dot(v_dyn_remain) > 0;
            
            std::vector<ContactPoint> grasp_mnps;
            std::vector<ContactPoint> mnp_fingertips;
            this->pw->world->getRobot()->getFingertipsOnObject(new_mnp_config, x, &mnp_fingertips);
            this->pw->world->getRobot()->Fingertips2PointContacts(mnp_fingertips, &grasp_mnps);
            
            Vector6d v_dyn = this->pw->QPVelocity_Quasidynamic(v, grasp_mnps, envs,
                    env_mode, f_o, T->translation_weight, T->angle_weight,1);

            bool isquasidynamic = v.dot(v_dyn) > 0.2*v.dot(v);
            
            bool ifCollide = this->pw->world->isRobotCollide(new_mnp_config); 

            if ((isremainbalance || isremaindynamic) && isquasidynamic && !ifCollide){
                new_config = new_mnp_config;
                return true;
            } 
        }
    // }

    return false;
}

void HPlanner::Extend_ObjectPath(int near_idx, Vector7d x_rand_){


    // extend all cs modes, best ss modes
    
    // steer goal
    Vector7d x_rand = steer_config(T->nodes[near_idx].config, x_rand_, this->eps_trans, this->eps_angle);
    // Vector7d x_rand = x_rand_;

    Vector6d v_star = compute_rbvel_body(T->nodes[near_idx].config, x_rand);

    Vector6d f_o = weight_w2o(T->nodes[near_idx].config, this->f_w);

    // contact mode enumeration
    
    MatrixXd A;
    VectorXd b;
    MatrixXd D;
    VectorXd d;
    this->pw->cons->NormalVelocityConstraints(T->nodes[near_idx].envs, &A, &b);
    this->pw->cons->TangentVelocityConstraints(T->nodes[near_idx].envs, &D, &d);

    if (T->nodes[near_idx].modes.size() == 0){
        T->nodes[near_idx].envs.clear();
        this->pw->world->getObjectContacts(&(T->nodes[near_idx].envs), T->nodes[near_idx].config);

        this->pw->cons->NormalVelocityConstraints(T->nodes[near_idx].envs, &A, &b);

        this->pw->cons->TangentVelocityConstraints(T->nodes[near_idx].envs, &D, &d);

        if (T->nodes[near_idx].envs.size() == 0){
            VectorXi m(0);
            T->nodes[near_idx].modes.push_back(m);
        } else{
            cs_mode_enumeration(A, &(T->nodes[near_idx].modes));
        }

    }

    // for every mode do forward integration
    Vector6d v_zero = Vector6d::Zero();
    double d_zero = dist_vel(v_zero, v_star, T->translation_weight, T->angle_weight);
    for(const auto& cs_mode: T->nodes[near_idx].modes){

        std::cout << "cs mode " << cs_mode.transpose() << std::endl;

        std::vector<VectorXi> modes;

        ss_mode_enumeration(A, D, cs_mode, &modes);

        ///
        // choose sliding modes
        double dv_best = d_zero;

        VectorXi mode_best;

        std::vector<VectorXi> mode_to_extend;
        {
            VectorXi all_sticking_mode(3*cs_mode.size());
            all_sticking_mode.setZero();
            all_sticking_mode.block(0, 0, cs_mode.size(),1) = cs_mode;
            Vector6d v = this->pw->EnvironmentConstrainedVelocity(v_star, T->nodes[near_idx].envs, all_sticking_mode);
            if (dist_vel(v, v_star, T->translation_weight, T->angle_weight) < d_zero){
                mode_to_extend.push_back(all_sticking_mode);
            }

        }
        for(const auto& mode: modes){
            Vector6d v = this->pw->EnvironmentConstrainedVelocity(v_star, T->nodes[near_idx].envs, mode);
            if (dist_vel(v, v_star, T->translation_weight, T->angle_weight) < dv_best)
            {
                dv_best = dist_vel(v, v_star, T->translation_weight, T->angle_weight);
                mode_best = mode;
            }
        }
        std::cout << "best mode " << mode_best.transpose() << std::endl; 
        
        if(dv_best < d_zero - 1e-4){
            mode_to_extend.push_back(mode_best);
        }

        /// choose sliding mode end

        for (const auto& mode: mode_to_extend){

            std::cout << "Extend mode: " << mode.transpose() << std::endl;
            
            std::vector<Vector7d> path;
            
            this->pw->ForwardIntegration_ObjectPath(
                                T->nodes[near_idx].config, x_rand, T->nodes[near_idx].envs, mode, 
                                T->translation_weight, T->angle_weight, &path);
        

            // if integration is successful
            if (path.size()>1){
                
                HNode new_node(path.back(), T->nodes[near_idx].dynamic_type);
                this->pw->world->getObjectContacts(&(new_node.envs), new_node.config);

                HEdge new_edge(mode, ifmanipulatorcollide, path);

                T->add_node(&new_node, near_idx, &new_edge);

            }
        }
        
    }
    return;

}

void HPlanner::Search(const RRTPlannerOptions& opts, Vector7d x_goal, double goal_thr, std::vector<int>* node_path, bool& ifsuccess, double& time_cpu_s, bool ifquasidynamic){
    
    this->options = new RRTPlannerOptions(opts);

    set_rand_seed(); //set random seed by cuurent time
    int goal_idx = -1;
    double goal_biased_prob = this->options->goal_biased_prob;

    std::clock_t c_start = std::clock();

    for(int kk = 0; kk <this->options->max_samples; kk++){
        std::cout << "iter: " << kk << std::endl;
        
        
        // bias sample toward the goal
        Vector7d x_rand;
        int near_idx;
        if (randd() > goal_biased_prob){
            Vector3d p_rand;
            Quaterniond q_rand;
            p_rand = sample_position(this->X_ub, this->X_lb);

            if (this->options->sampleSO3){
                q_rand = generate_unit_quaternion();
            } else {
                q_rand = sample_rotation(this->options->rotation_sample_axis);
            }

            x_rand << p_rand[0], p_rand[1], p_rand[2],
                      q_rand.x(), q_rand.y(), q_rand.z(), q_rand.w();
            near_idx = T->nearest_neighbor(x_rand);
        } else {
            x_rand = x_goal;
            // near_idx = T->nearest_neighbor(x_rand);
            near_idx = T->nearest_unextended_to_goal(x_goal);
            if (near_idx != 0){
                T->nodes[near_idx].is_extended_to_goal = true;
            }
            
        }
        if (kk<1){
            x_rand = x_goal;
            near_idx = 0;
        }
        
        std::cout << "x_rand:" << x_rand.transpose() << ", near idx: "<< near_idx << std::endl;

        this->Extend_ObjectPath(near_idx, x_rand);

        int goal_near_idx = T->nearest_neighbor(x_goal);
        double dd = T->dist(T->nodes[goal_near_idx].config, x_goal); 

        // continue to find object path if the goal is not reached yet
        if (dd > goal_thr)
        {
            continue;
        }
        
        printf("Found goal node in %d samples. \n", kk + 1);
        std::cout << T->nodes[goal_near_idx].manipulator_config.transpose() << std::endl;
        goal_idx = goal_near_idx;
        
        std::vector<int> current_node_path;
        T->backtrack(goal_idx, current_node_path);
        std::reverse(current_node_path->begin(), current_node_path->end());

        bool if_subplan_success = false;
        int max_n_feasible = 0;

        ///* TODO: subplanner here : plan manipulator contacts and check if the plan is feasible *///

        for (int k_subplan = 0; k_subplan < 15; k_subplan++){
            int n_feasible = 0;
            std::vector<Vector7d> object_configs; 
            std::vector<VectorXd> mnp_configs;
            if_subplan_success = this->GenerateFingerPlan(current_node_path, n_feasible,
                                        &object_configs, &mnp_configs);
            if (if_subplan_success){
                ///* TODO: show the plan here *///
                break;
            }
            else {
                if (max_n_feasible < n_feasible){
                    max_n_feasible = n_feasible;
                }
            }
        }

        if (if_subplan_success){
            break;
        }
        else {
            // kill infeasible nodes
            goal_idx = -1;

            current_node_path.erase(current_node_path.begin(), current_node_path.begin() + max_n_feasible);

            for (int kill_idx:current_node_path){
                this->T->nodes[kill_idx].is_alive = false;
            }
            // Debug: this block shoudn't be needed
            for (int kill_idx = 1; kill_idx < this->T->nodes.size(); kill_idx++){

                if(!this->T->nodes[this->T->nodes[kill_idx].parent].is_alive){
                    this->T->nodes[kill_idx].is_alive = false;
                }
            }

            printf("Cannot find a finger plan. Killed infeasible nodes. Max length of finger plan: %d.\n", max_n_feasible);

        }
        

    }

    // time
    std::clock_t c_end = std::clock();

    if (goal_idx!=-1){
        ifsuccess = true;
        printf("GOAL REACHED! \n");
    } else {
        ifsuccess = false;
        goal_idx =  T->nearest_neighbor(x_goal);
        std::cout << "GOAL NOT REACHED. Dist: " << T->dist(T->nodes[goal_idx].config, x_goal)
        //  << ". Closest config: " << T->nodes[goal_idx].config.transpose() 
         << std::endl;
        
    }
    
    time_cpu_s = double(c_end-c_start) / double(CLOCKS_PER_SEC);
    std::cout << "CPU time used: " << time_cpu_s << " s\n";
    std::cout << "Goal idx " << goal_idx << std::endl;
    // backtrack
    T->backtrack(goal_idx, node_path);
    std::reverse(node_path->begin(), node_path->end());

    return;
}

void HPlanner::VisualizePath(const std::vector<int>& node_path){
    // visualization
    std::vector<Vector7d> object_poses;
    std::vector<VectorXd> mnp_configs;

    for(auto& k:node_path){
        if(k==0){continue;}
        std::vector<Vector7d> path = T->edges[T->nodes[k].edge].path;
        VectorXd mnp_config = T->nodes[k].manipulator_config;
        // std::cout << mnp_config.transpose() << std::endl;
        for (auto x: path){
            object_poses.push_back(x);
            mnp_configs.push_back(mnp_config);
        }
    }

    this->pw->world->setPlaybackTrajectory(object_poses, mnp_configs);

}

bool HPlanner::GenerateFingerPlan(const std::vector<int> node_path, int& feasible_nodes, 
    std::vector<Vector7d>* object_configs, std::vector<VectorXd>* mnp_configs)
{

    // TODO: do we need to pass dynamic types to the tree nodes?
    std::vector<int> dynamic_types;
    // the initial node is always quasistatic
    dynamic_types.push_back(QUASISTATIC);

    VectorXd mnp_config = T->nodes[node_path[0]].manipulator_config;

    for (feasible_nodes = 1; feasible_nodes < node_path.size(); feasible_nodes++){
        
        int start_node_idx = node_path[feasible_nodes-1];
        int goal_node_idx = node_path[feasible_nodes];

        // int dyn_type = T->nodes[goal_node_idx].dynamic_type;
        int dyn_type = (dynamic_types.back() == DYNAMIC) ? DYNAMIC : QUASISTATIC;
                
        VectorXi mode = T->edges[T->nodes[goal_node_idx].edge].mode;
        
        std::vector<ContactPoint> envs_pre = 
                                T->nodes[T->nodes[goal_node_idx].parent].envs;
            
        std::vector<Vector7d> node_obj_path;
        std::vector<VectorXd> node_mnp_path;

        // Forward check
        // add new nodes while doing forward check

        // TODO: forward check for all three types

        bool is_forwarded = false;

        if (dyn_type == QUASISTATIC){
            // TODO
            node_obj_path = T->edges[T->nodes[goal_node_idx].edge].path;
            is_forwarded = ForwardCheck_Quasistatic(T->nodes[start_node_idx].config,T->nodes[goal_node_idx].config, 
                mnp_config, envs_pre, mode, this->f_w, &node_obj_path, &node_mnp_path);
            
            if(!is_forwarded){
                dyn_type = QUASIDYNAMIC;
            }
        }

        if (dyn_type == QUASIDYNAMIC){

            // TODO
            is_forwarded = ForwardCheck_Quasidynamic(T->nodes[start_node_idx].config,T->nodes[goal_node_idx].config, 
                mnp_config, envs_pre, mode, this->f_w, &node_obj_path, &node_mnp_path);

            if(!is_forwarded){
                dyn_type = DYNAMIC;
            }

        }

        if (dyn_type == DYNAMIC){

            // TODO
            is_forwarded = ForwardCheck_Dynamic(T->nodes[start_node_idx].config,T->nodes[goal_node_idx].config, 
                mnp_config, envs_pre, mode, this->f_w, &node_obj_path, &node_mnp_path);

        } 

        
        object_configs->insert(object_configs->end(), node_obj_path.begin(), node_obj_path.end());
        mnp_configs->insert(mnp_configs->end(), node_mnp_path.begin(), node_mnp_path.end());
        dynamic_types.push_back(dyn_type);

        if(!is_forwarded){
            return false;
        }

        mnp_config = mnp_configs->back();
    }

    return true;
}

bool HPlanner::ForwardCheck_Quasistatic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path)
{
    if (mnp_path->size() > 0){
        mnp_path->clear();
    }
    VectorXd mnp_config = pre_mnp_config;

    int idx = 0;
    for (idx = 0; idx < path->size(); idx++){
        
        Vector7d x = path->at(idx);

        Vector6d f_o = weight_w2o(x, this->f_w);

        // collision detection
        std::vector<ContactPoint> envs;
        this->pw->world->getObjectContacts(&envs, x);

        if(envs.size() != env_mode_.size()/3){
            printf("Error in HPlanner::ForwardCheck_Quasistatic(). Collision detection inconsistent with prespecified mode.\n");            
            break;
        }

        
        bool if_relocate = (randd() < 0.1) ? true: false;
        bool is_balance = false;
        
        if (mnp_config.size() == 0){
            if_relocate = true;
        } else {

            std::vector<ContactPoint> grasp_mnps;
            std::vector<ContactPoint> mnp_fingertips;
            this->pw->world->getRobot()->getFingertipsOnObject(
                mnp_config, x, &mnp_fingertips);
            this->pw->world->getRobot()->Fingertips2PointContacts(mnp_fingertips, &grasp_mnps);

            is_balance = this->pw->isForceBalance(grasp_mnps, envs_, 
                                        env_mode_, f_o);
            
            if (!is_balance){
                if_relocate = true;
            }
        }

        if (if_relocate){
            bool ifchanged;

            VectorXd new_mnp_config;
                
            ifchanged = this->randomRelocateFingers(mnp_config, x, 
                            envs, env_mode_, f_o, new_mnp_config);
        
            if(!ifchanged && !is_balance){ 
                // cannot find good mnp locations, continue to the next mode
                printf("Cannot find good mnp locations\n");
                break;
            } else {
                mnp_config = new_mnp_config;
            }
        }

        mnp_path->push_back(mnp_config);
    }

    if (mnp_path->size() == path->size()){
        return true;
    }
    path->erase(path->begin()+mnp_path->size(), path->end());
    return false;
}

bool HPlanner::ForwardCheck_Quasidynamic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path)
{
    if (mnp_path->size() > 0){
        mnp_path->clear();
    }

    if (path->size() > 0){
        path->clear();
    }

    double thr = 1e-3;
    double h = 0.04;
    int max_counter = 200;

    bool ifmanipulatorcollide = false;

    Vector7d x = x_start;
    this->world->updateObjectPose(x);

    VectorXi env_mode = env_mode_;

    VectorXd mnp_config = pre_mnp_config;

    path->push_back(x);
    mnp_path->push_back(mnp_config);

    std::vector<ContactPoint> envs;
    std::vector<ContactPoint> mnps;
    std::vector<ContactPoint> fingertips;
    this->world->getRobot()->getFingertipsOnObject(mnp_config, x, &fingertips);
    this->world->getRobot()->Fingertips2PointContacts(fingertips, &mnps);
    envs = envs_;

    std::vector<ContactPoint> envs_pre;
    envs_pre = envs;

    Vector6d v_b_pre = Vector6d::Zero();

    int counter;
    int delete_c=0;

    bool ifChangeFingertips = false;


    for (counter=0; counter < max_counter; counter++) 
    {
        Vector6d v_star = compute_rbvel_body(x, x_goal);

        if (v_star.norm() < thr){
            std::cout << "v_star < thr : " << v_star.transpose() << std::endl;
            break;
        }

        Matrix4d T = pose2SE3(x);
        Matrix6d Adg = SE32Adj(T);

        Matrix4d T_;
        T_.setIdentity();
        T_.block(0,0,3,3) = T.block(0,0,3,3);
        Vector6d f_ext_o = SE32Adj(T_).transpose()*f_ext_w;

        Vector6d v_b;

        // change fingertips if needed: begin
        bool isforcefeasible = false;

        if(mnp_config.size()==0){
            // if finger not assigned
            ifChangeFingertips = true;
        }
        else {
            bool ifmnpcollide = this->pw->world->isRobotCollide(mnp_config);
            
            if (ifmnpcollide){
                ifChangeFingertips = true;
            } else {
                // check quasidynamic feasibility

                // bool isforcefeasible = this->pw->isForceBalance(fingertips, envs, mode, f_o);
                Vector6d v_dyn = this->pw->QPVelocity_Quasidynamic(v_star, mnps, envs,
                env_mode, f_ext_o, this->T->translation_weight, this->T->angle_weight,1);

                isforcefeasible = v_star.dot(v_dyn) > 0.2*v_star.dot(v_star);

                if (!isforcefeasible){
                    // if yes no need to change fingertips
                    if (v_star.norm() < 0.1){
                        return true;
                    } else {
                        ifChangeFingertips = true;
                    }
                    
                } else {
                    v_b = v_dyn;
                }
            }
        }

        if (randd() < 0.1){
            ifChangeFingertips = true;
        }

        if (ifChangeFingertips){

            VectorXd new_mnp_config;

            bool ifchanged = this->randomRelocateFingers_Quasidynamic(mnp_config, x, v_star,
                            envs, env_mode, f_ext_o, new_mnp_config);
        
            if(!ifchanged){
                if (!isforcefeasible){
                // cannot find good mnp locations, continue to the next mode
                    printf("ForwardCheck_Quasidynamic: Cannot find good mnp locations\n");
                    break;
                }
            } else {
                mnp_config = new_mnp_config;
                mnps.clear();
                this->pw->manipulator->pointContactApproximation(mnp_config, &mnps);
            }

            ifChangeFingertips = false;
        } 
        // Change Fingertips If Needed: End

        if (v_b.norm() < 1e-3){
            v_b = this->pw->QPVelocity_Quasidynamic(v_star, mnps, envs, env_mode, f_ext_o, this->T->translation_weight, this->T->angle_weight, 1);
        }
        
        if (v_b.norm() < thr ){
            std::cout << "v_b < thr : " << v_b.transpose() << std::endl;
            break;
        }

        if ((v_b_pre.transpose()*v_b)[0] < -1e-5){
            printf("v_b back and forth. \n");
            if (v_star.norm() < 0.1){
                return true;
            } else {
                ifChangeFingertips = true;
            }

        }

        steer_velocity(v_b, h, this->charac_len);

        // integrate v
        Vector7d x_new = SE32pose(T*se32SE3(v_b));

        // check manipulator collision (break may happen)
        this->world->updateObjectPose(x_new);

        // check penetration & interpolation (break may happen)
        envs.clear();
        this->world->getObjectContacts(&envs, x_new);

        // velocity correction
        if (envs.size()!=0 && envs.size() == env_mode.size()/3 && (ifNeedVelocityCorrection(env_mode, envs))){
            // std::cout << "velocity correction "<< counter << std::endl;
            Vector6d v_corr = VelocityCorrection(envs);
            x_new = SE32pose(pose2SE3(x_new)*se32SE3(v_corr));
            envs.clear();
            this->world->getObjectContacts(&envs, x_new);
        }

        if (envs.size()>env_mode.size()/3){
            Vector6d v_corr = VelocityCorrection(envs);
            x_new = SE32pose(pose2SE3(x_new)*se32SE3(v_corr));;
            path->push_back(x_new);
            x = x_new;
            envs.clear();
            this->world->getObjectContacts(&envs, x_new);
            printf("Made new contacts! \n");
            print_contacts(envs);
            break;
        }

        // update contact mode if needed
        if (envs.size()<env_mode.size()/3){
            VectorXi remain_idx = track_contacts_remain(envs_pre, envs);
            if (envs.size()!=0 && remain_idx.size() == 0){
                printf("contact track fails.\n");
                break;
            }
            if (ifContactingModeDeleted(env_mode, remain_idx)){
                if (h < 0.004/3){
                    // TODO: need to fix this
                    delete_c++;
                    if (delete_c > 4)
                        break;
                    printf("Collision Detection delelte contacting mode \n");
                    Vector6d v_corr = recoverContactingContacts(envs_pre, env_mode, remain_idx);
                    x_new = SE32pose(pose2SE3(x_new)*se32SE3(v_corr));
                    envs.clear();
                    this->world->getObjectContacts(&envs, x_new);
                    if(envs.size() <= env_mode.size()/3)
                        continue;
                    else
                        break;
                }
                h = h/1.5;
                envs = envs_pre;
                continue;
            } else {
                env_mode = deleteModebyRemainIndex(env_mode, remain_idx);
            } 
        }

        x = x_new;
        envs_pre = envs;
        v_b_pre = v_b;
        
        path->push_back(x);
        mnp_path->push_back(mnp_config);

        if(counter == max_counter - 1){
            printf("Reach the end.\n");
        }

    }

    double d_goal = (compute_rbvel_body(x, x_goal)).norm();
    return d_goal < 0.1;
}

bool HPlanner::ForwardCheck_Dynamic(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& pre_mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    std::vector<Vector7d>* path, std::vector<VectorXd>* mnp_path)
{
    if (mnp_path->size() > 0){
        mnp_path->clear();
    }

    if (path->size() > 0){
        path->clear();
    }

    return false;
}