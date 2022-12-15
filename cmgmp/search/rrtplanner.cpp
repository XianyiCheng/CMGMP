#include "rrtplanner.h"

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

Vector7d steer_config_directional(Vector7d x_near, Vector7d x_rand, Vector6d v, double epsilon_translation, double epsilon_angle){
    
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


RRTPlanner::RRTPlanner(PlanningWorld* pw_){
    this->pw = pw_;
}

void RRTPlanner::Initialize(Vector3d q_lb, Vector3d q_ub, Vector7d q_start, Vector6d f_w, std::vector<ContactPoint>& surface, 
    double angle_weight, double translation_weight, double eps_trans, double eps_angle)
{
    // angle_weight and translation_weight: heuristics of how important the dimensions are. 

    // this->T = new Tree(angle_weight, translation_weight);
    this->T = new Tree(angle_weight, this->pw->charac_len*translation_weight);
    this->X_ub = q_ub;
    this->X_lb = q_lb;
    this->f_w = f_w;
    this->eps_trans = eps_trans;
    this->eps_angle = eps_angle;

    this->object_surface = surface;

    VectorXd empty_mnps(0);
    Node start_node(q_start, empty_mnps);
    this->pw->world->getObjectContacts(&(start_node.envs), start_node.config);
    T->initial_node(&start_node); 

}

void RRTPlanner::SetInitialManpulatorConfig(const VectorXd& config){
    T->nodes[0].manipulator_config = config;
}

bool RRTPlanner::randomRelocateFingers(const VectorXd& mnp_config, Vector7d x,
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
    VectorXd& new_config)
{
    // std::vector<ContactPoint> fingertips;
    // this->pw->world->getRobot()->getFingertips(mnp_config, &fingertips);

    this->pw->world->updateObjectPose(x);

    VectorXd new_mnp_config(0);
    // TODO: this can be again formulated as a search problem
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

bool RRTPlanner::randomRelocateFingers_Quasidynamic(const VectorXd& mnp_config, Vector7d x, Vector6d v,
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode,  Vector6d f_o, 
    VectorXd& new_config)
{
    // std::vector<ContactPoint> fingertips;
    // this->pw->world->getRobot()->getFingertips(mnp_config, &fingertips);
    this->pw->world->updateObjectPose(x);
    VectorXd new_mnp_config(0);
    // TODO: this can be again formulated as a search problem
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
            
            bool ifCollide = this->pw->world->isRobotCollide(new_mnp_config); // TODO

            if ((isremainbalance || isremaindynamic) && isquasidynamic && !ifCollide){
                new_config = new_mnp_config;
                return true;
            } 
        }
    // }

    return false;
}

void RRTPlanner::Extend_QPCC(int near_idx, Vector7d x_rand_){


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
        // std::vector<ContactPoint> envs;
        // envs.push_back(T->nodes[near_idx].envs[2]);
        // envs.push_back(T->nodes[near_idx].envs[0]);
        // envs.push_back(T->nodes[near_idx].envs[3]);
        // envs.push_back(T->nodes[near_idx].envs[1]);

        this->pw->cons->NormalVelocityConstraints(T->nodes[near_idx].envs, &A, &b);
        // this->pw->cons->NormalVelocityConstraints(envs, &A, &b);
        this->pw->cons->TangentVelocityConstraints(T->nodes[near_idx].envs, &D, &d);
        // this->pw->cons->TangentVelocityConstraints(envs, &D, &d);
        // std::cout << "A\n" << A <<std::endl;
        // std::cout << "D\n" << D <<std::endl;

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

        // print_contacts(T->nodes[near_idx].envs);
        
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
            // std::cout << "mode: " << all_sticking_mode.transpose() <<  " velocity: " << v.transpose() << std::endl;
            if (dist_vel(v, v_star, T->translation_weight, T->angle_weight) < d_zero){
                mode_to_extend.push_back(all_sticking_mode);
            }

        }
        for(const auto& mode: modes){
            Vector6d v = this->pw->EnvironmentConstrainedVelocity(v_star, T->nodes[near_idx].envs, mode);
            // std::cout << "mode: " << mode.transpose() <<  " velocity: " << v.transpose() << std::endl;
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
            // check if we can try direct forward integration
            bool ifChangeFingertips = false;

            VectorXd mnp_config;

            if(T->nodes[near_idx].manipulator_config.size()==0){
                ifChangeFingertips = true;
            } else {
                bool ifmnpcollide;
                if(T->nodes[near_idx].parent == -1){
                    ifmnpcollide = false;
                } else {
                    ifmnpcollide = T->edges[T->nodes[near_idx].edge].manipulator_collide;
                }
                if (ifmnpcollide){
                    ifChangeFingertips = true;
                } else {
                    // check if there is force balance
                    std::vector<ContactPoint> grasp_mnps;
                    std::vector<ContactPoint> mnp_fingertips;
                    this->pw->world->getRobot()->getFingertipsOnObject(
                        T->nodes[near_idx].manipulator_config, T->nodes[near_idx].config, &mnp_fingertips);
                    this->pw->world->getRobot()->Fingertips2PointContacts(mnp_fingertips, &grasp_mnps);

                    bool isbalance = this->pw->isForceBalance(grasp_mnps, T->nodes[near_idx].envs, 
                                    mode, f_o);
                    
                    if (!isbalance){
                        // if yes no need to change fingertips
                        ifChangeFingertips = true;
                    }
                }
            }

            // if no try change manipulator contacts and then do forward integration
            if (ifChangeFingertips){
                bool ifchanged;
                 
                ifchanged = this->randomRelocateFingers(T->nodes[near_idx].manipulator_config, T->nodes[near_idx].config, 
                                T->nodes[near_idx].envs, mode, f_o, mnp_config);
            
                if(!ifchanged){ 
                    // cannot find good mnp locations, continue to the next mode
                    printf("Cannot find good mnp locations\n");
                    continue;
                } 
            } else {
                mnp_config = T->nodes[near_idx].manipulator_config;
            } 

            std::vector<Vector7d> path;
            bool ifmanipulatorcollide; 
            
            ifmanipulatorcollide = this->pw->ForwardIntegration(
                                T->nodes[near_idx].config, x_rand, mnp_config, T->nodes[near_idx].envs, mode, this->f_w, 
                                T->translation_weight, T->angle_weight, &path);
        

            // if integration is successful
            if ((path.size()>1 && !ifmanipulatorcollide) || (path.size() > 2)){
                
                Node new_node(path.back(), mnp_config);
                this->pw->world->getObjectContacts(&(new_node.envs), new_node.config);

                Edge new_edge(mode, ifmanipulatorcollide, path);

                T->add_node(&new_node, near_idx, &new_edge);

            }
        }
        
    }
    return;

}

void RRTPlanner::Extend_QPCC_Quasidynamic(int near_idx, Vector7d x_rand_){

    bool ifquasidynamic = true;

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
            // std::cout << "mode: " << all_sticking_mode.transpose() <<  " velocity: " << v.transpose() << std::endl;
            if (dist_vel(v, v_star, T->translation_weight, T->angle_weight) < d_zero){
                mode_to_extend.push_back(all_sticking_mode);
            }

        }
        for(const auto& mode: modes){
            Vector6d v = this->pw->EnvironmentConstrainedVelocity(v_star, T->nodes[near_idx].envs, mode);
            // std::cout << "mode: " << mode.transpose() <<  " velocity: " << v.transpose() << std::endl;
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

            Vector6d v_mode = this->pw->EnvironmentConstrainedVelocity(v_star, T->nodes[near_idx].envs, mode);
            // check if we can try direct forward integration
            bool ifChangeFingertips = false;

            VectorXd mnp_config;

            // if(T->nodes[near_idx].manipulator_config.size()==0){
            //     ifChangeFingertips = true;
            // } else 
            {
                bool ifmnpcollide;
                if(T->nodes[near_idx].parent == -1){
                    ifmnpcollide = false;
                } else {
                    ifmnpcollide = T->edges[T->nodes[near_idx].edge].manipulator_collide;
                }
                if (ifmnpcollide){
                    ifChangeFingertips = true;
                } else {
                    // check if there is force balance
                    std::vector<ContactPoint> grasp_mnps;
                    std::vector<ContactPoint> mnp_fingertips;
                    this->pw->world->getRobot()->getFingertipsOnObject(
                        T->nodes[near_idx].manipulator_config, T->nodes[near_idx].config, &mnp_fingertips);
                    this->pw->world->getRobot()->Fingertips2PointContacts(mnp_fingertips, &grasp_mnps);

                    Vector6d v_dyn = this->pw->QPVelocity_Quasidynamic(v_mode, grasp_mnps, T->nodes[near_idx].envs,
                        mode, f_o, T->translation_weight, T->angle_weight,1);
                    
                    // Vector6d v_static = this->pw->QPVelocity(v_star, fingertips, T->nodes[near_idx].envs,
                    //     mode, f_o, T->translation_weight, T->angle_weight);
                    
                    // bool isforcebalance = this->pw->isForceBalance(fingertips, T->nodes[near_idx].envs,
                    //     mode, f_o);

                    bool isdynamicfeasible;
                    if (v_dyn.norm() < 1e-4){
                        isdynamicfeasible = false;
                    } else {
                        isdynamicfeasible = v_mode.dot(v_dyn*(1/v_mode.norm())) > 0.2*v_mode.dot(v_mode*(1/v_mode.norm()));

                    }
                    
                    if (!isdynamicfeasible){
                        // if yes no need to change fingertips
                        ifChangeFingertips = true;
                    }
                }
            }

            // if no try change manipulator contacts and then do forward integration
            if (ifChangeFingertips){
                bool ifchanged;
                 
                ifchanged = this->randomRelocateFingers_Quasidynamic(T->nodes[near_idx].manipulator_config, T->nodes[near_idx].config, v_mode, 
                                T->nodes[near_idx].envs, mode, f_o, mnp_config);
            
                if(!ifchanged){ 
                    // cannot find good mnp locations, continue to the next mode
                    printf("Cannot find good mnp locations\n");
                    continue;
                } 
            } else {
                mnp_config = T->nodes[near_idx].manipulator_config;
            } 

            std::vector<Vector7d> path;
            bool ifmanipulatorcollide; 
            
            ifmanipulatorcollide = this->pw->ForwardIntegration(
                                T->nodes[near_idx].config, x_rand, mnp_config, T->nodes[near_idx].envs, mode, this->f_w, 
                                T->translation_weight, T->angle_weight, &path, ifquasidynamic);
        

            // if integration is successful
            if (path.size()>1){
                
                Node new_node(path.back(), mnp_config);
                this->pw->world->getObjectContacts(&(new_node.envs), new_node.config);

                Edge new_edge(mode, ifmanipulatorcollide, path);

                T->add_node(&new_node, near_idx, &new_edge);

            }
        }
        
    }
    return;

}

void RRTPlanner::Search(const RRTPlannerOptions& opts, Vector7d x_goal, double goal_thr, std::vector<int>* node_path, bool& ifsuccess, double& time_cpu_s, bool ifquasidynamic){
    
    this->options = new RRTPlannerOptions(opts);

    set_rand_seed(); //set random seed by cuurent time
    int goal_idx = -1;
    // double goal_thr = 3.14*1/180; // TODO: adjust threshold
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

        if (ifquasidynamic) {
            this->Extend_QPCC_Quasidynamic(near_idx, x_rand);
        } else {
            this->Extend_QPCC(near_idx, x_rand);
        }

        
        /*
        // extend function
        switch(this->options->method) {
            case METHOD_QPCC:
                this->Extend_QPCC(near_idx, x_rand);
                break; 
            case METHOD_LCP:
                this->Extend_LCP(near_idx, x_rand);
                break; 
            case METHOD_ALL:
                this->Extend_All(near_idx, x_rand);
                break; 
            case METHOD_NONE:
                this->Extend_None(near_idx, x_rand);
                break; 
            case METHOD_QPCC_BACKWARD:
                this->Extend_QPCC_BACKWARD(near_idx, x_rand);
                break; 
            default : 
                this->Extend_QPCC(near_idx, x_rand);
        }
        */

        int goal_near_idx = T->nearest_neighbor(x_goal);
        double dd = T->dist(T->nodes[goal_near_idx].config, x_goal); 
        if (dd <= goal_thr)
        {    
            printf("Found goal node in %d samples. \n", kk + 1);
            std::cout << T->nodes[goal_near_idx].manipulator_config.transpose() << std::endl;
            goal_idx = goal_near_idx;
            break;
        }

        if (kk % 100 == 0) {
            std::cout << "Iteration: " << kk << " Nodes Expanded: " << T->nodes.size() << std::endl;
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
         << ". Closest config: " << T->nodes[goal_idx].config.transpose() << std::endl;
        
    }

    // std::cout << "Nearest Mode: " <<  T->edges[T->nodes[goal_idx].edge].mode.transpose() << std::endl;
    
    time_cpu_s = double(c_end-c_start) / double(CLOCKS_PER_SEC);
    std::cout << "CPU time used: " << time_cpu_s << " s\n";
    std::cout << "Goal idx " << goal_idx << std::endl;
    // backtrack
    T->backtrack(goal_idx, node_path);
    std::reverse(node_path->begin(), node_path->end());

    return;
}

void RRTPlanner::VisualizePath(const std::vector<int>& node_path){
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