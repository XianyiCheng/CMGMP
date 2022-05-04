#include "PlanningWorld.h"

#ifndef SAMPLE_H
#define SAMPLE_H
    #include "../utilities/sample.h"
#endif

#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "../contacts/contact_kinematics.h"
#endif 

#include "../utilities/eiquadprog.hpp"
#include <math.h>


#define OBJECT_USER_ID 1

// btTransform pose2btTransform(const Vector7d& x){
//     btTransform tr(btQuaternion(x(3), x(4), x(5), x(6)), btVector3(x(0), x(1), x(2)));
//     return tr;
// }

bool ifConstraintsSatisfied(const VectorXd& x, const MatrixXd A, const VectorXd& b, const MatrixXd G, const VectorXd& h){

    double Ax_b = (A*x - b).cwiseAbs().sum();
    if (Ax_b > 1e-3){ 
        return false;
    }

    VectorXd g = G*x - h;
    for(int i = 0; i < g.size(); i++){
        if (g[i] < -1e-3){
            return false;
        }
    }

    return true;
}

void steer_velocity(Vector6d& x, double h, double cl = 1.0){

    Vector6d xx = x;
    xx.head(3) = xx.head(3)*cl;
    if (xx.norm() > h){
        xx = (h/xx.norm())*xx;
        double xx_norm = xx.tail(3).norm();
        double x_norm = x.tail(3).norm();
        if (xx_norm > 1e-4 && x_norm > 1e-4){
            x = x*(xx.tail(3).norm()/x.tail(3).norm());
        } else {
            x = xx/cl;
        }
    }
    return;
}



bool ifCollide(const std::vector<ContactPoint>& pts){
    double thr = -0.04;
    for (const auto& pt: pts){
        if (pt.d < thr){
            return true;
        }
    }
    return false;
}

VectorXi deleteSeparatingMode(const VectorXi& mode){
    
    VectorXi m(mode.size());
    int n_cpts = 0;
    int n_pts = mode.size()/3;

    for (int i = 0; i < n_pts; i++){
        int cs_mode = mode(i);
        if (cs_mode == 0){ // contacting
            n_cpts+=1;
        }
    }
    int k = 0;
    for (int i = 0; i < n_pts; i++){
        int cs_mode = mode(i);
        if (cs_mode == 0){ // contacting
            m(k) = cs_mode;
            m.block(n_cpts+2*k,0,2,1) = mode.block(n_pts+2*i,0,2,1);
            k += 1;
        }
    }
    m.conservativeResize(3*n_cpts);
    return m;

}

VectorXi deleteModebyRemainIndex(const VectorXi& mode, const VectorXi& remain_idx){
    
    int n_cpts = remain_idx.size();
    int n_pts = int(mode.size()/3);
    VectorXi m(3*n_cpts);
    for (int i = 0; i < n_cpts; i++){
        m[i] = mode[remain_idx[i]];
        m.block(n_cpts+2*i,0,2,1) = mode.block(n_pts+2*remain_idx[i],0,2,1);
    }
    return m;
}

bool ifContactingModeDeleted(const VectorXi& mode, const VectorXi& remain_idx){
    for (int i = 0; i < int(mode.size()/3); i++){
        if(mode[i] == 0){
            bool ifremained = false;
            for (int k = 0; k < remain_idx.size(); k++){
                if (i == remain_idx[k]){
                    ifremained = true;
                    break;
                }
            }
            if(!ifremained){
                return true;
            }
        }
    }
    return false;
}

double CollisionInterpolation(const Vector6d& v, const std::vector<ContactPoint>& pts){
    
    double d_min = 0;
    Vector3d p;
    Vector3d n;

    for (const auto& pt: pts){
        if (abs(pt.d) > abs(d_min)){
            d_min = pt.d;
            p = pt.p;
            n = pt.n;
        }
    }
    Vector3d vel =  v.block(0,0,3,1);
    Vector3d omega = v.block(3,0,3,1);
    Vector3d v_p_max;
    v_p_max = omega.cross(p) + vel;
    double k = 0;
    double a = (v_p_max.transpose()*n)(0);
    if (std::abs(a) >= std::abs(d_min)){
        if (d_min > 0){
            k = 1 + (std::abs(d_min)-0.005) / std::abs(a);
        } else{
            k = 1 - (std::abs(d_min)-0.005) / std::abs(a);
        }
    }
    return k;
}

bool ifNeedVelocityCorrection(VectorXi mode, const std::vector<ContactPoint>& pts){
    double thr = 0.03;
    for (int i = 0; i < pts.size(); i++){
        
        if ((abs(pts[i].d) > thr) && mode[i] == 0){
            return true;
        }
    }
    return false;
}

Vector6d VelocityCorrection(const std::vector<ContactPoint>& pts){
    
    int n_pts = pts.size();
    double violation = -1e-4;

    Vector6d z_axis;
    z_axis << 0,0,1,0,0,0;
    MatrixXd N(n_pts, 6);
    VectorXd d(n_pts);

    for(int i = 0; i < n_pts; i++){
        Matrix6d Adgco = contact_jacobian(pts[i].p, pts[i].n);
        N.block(i,0,1,6) = z_axis.transpose()*Adgco;
        
        d[i] = -(pts[i].d - violation);
        
        
    }

    Matrix6d I;
    I.setIdentity();
    Matrix6d G;
    G = N.transpose()*N + 0.001*I;
    Vector6d g0;
    g0 = -2*(d.transpose()*N).transpose();

    // MatrixXd A(6,0);
    // VectorXd b(0);
    // MatrixXd C(6,0);
    // VectorXd e(0);

    Vector6d x;
    x = -G.inverse()*g0/2;

    // double f = solve_quadprog(G, g0, A, b,  C, e, x);

    // if (f > 0.0){
    //     x.setZero();
    // }


    return x;
}

bool contactTrack(ContactPoint pt0, ContactPoint pt1){
    if(((pt0.p - pt1.p).norm() < 0.1) && ((pt0.n.transpose()*pt1.n)[0] > 0.85)){
        return true;
    } else {
        // std::cout << "d: " << (pt0.p - pt1.p).norm() << " ";
        return false;
    }

}

VectorXi track_contacts_remain(const std::vector<ContactPoint>& pts, const std::vector<ContactPoint>& pts_new){
    VectorXi remain_idx(pts_new.size());
    int i = 0;
    int j = 0;
    while((i < pts.size()) && (j < pts_new.size())){
        if(contactTrack(pts[i], pts_new[j])){
            remain_idx[j] = i;
            j++;
        }
        i++; 
    }
    if (j < pts_new.size()){
        VectorXi empty_idx(0);
        return empty_idx;
    }
    return remain_idx;
}

void deleteExtraContacts(const std::vector<ContactPoint>& pts0, std::vector<ContactPoint>& pts){
    std::vector<ContactPoint> pts2;
    VectorXi track(pts0.size());
    track.setZero();

    for (auto& pt: pts){
        bool iftracked = false;
        for(int i = 0; i < pts0.size(); i++){
            if(track[i] == 1)
                continue;
            ContactPoint pt0 = pts0[i];
            iftracked = contactTrack(pt0, pt);
            if (iftracked){
                track[i] = 1;
                break;
            }
                
        }
        if (iftracked)
            pts2.push_back(pt);
    }
    pts = pts2;
}

bool simplify_line_contacts(const std::vector<ContactPoint>& pts, std::vector<ContactPoint>* pts_update){
    
    if (pts.size() <= 2){
        return false;
    }

    double thr = 5e-2;
    Vector3d vec = pts[0].p - pts[1].p;
    vec = vec/vec.norm();
    int idx = 1;
    double d = vec.norm();
    // check if pts are in the same line
    for (int i = 2; i < pts.size(); i++){
        Vector3d vec1 = pts[i].p - pts[0].p;
        vec1 = vec1/vec1.norm();
        double cross_product = vec.cross(vec1).norm();
        if (cross_product > thr){
            // not in the same line 
            return false;
        } else {
            double dd = (pts[i].p - pts[0].p).norm();
            if (dd > d){
                d = dd;
                idx = i;
            }
        }
    }
    pts_update->push_back(pts[0]);
    pts_update->push_back(pts[idx]);
    return true;
}

bool same_line_update(const std::vector<ContactPoint>& pts, 
const std::vector<ContactPoint>& pts_new, std::vector<ContactPoint>* pts_update)
{
    
    if (pts.size() != 2){
        return false;
    }

    double thr = 1e-3;
    Vector3d vec = pts[0].p - pts[1].p;
    double d = vec.norm();
    ContactPoint pt = pts[1];
    // check if pts_new are in the same line
    for (int i = 0; i < pts_new.size(); i++){
        if ((vec.cross(pts_new[i].p - pts[0].p)).norm() > thr){
            // not in the same line 
            return false;
        } else {
            double dd = (pts_new[i].p - pts[0].p).norm();
            if (dd > d){
                d = dd;
                pt = pts_new[i];
            }
        }
    }
    pts_update->push_back(pts[0]);
    pts_update->push_back(pt);
    return true;

}

Vector6d recoverContactingContacts(const std::vector<ContactPoint>& pts, const VectorXi& mode, const VectorXi& remain_idx){
   std::vector<ContactPoint> envs;
   for (int i = 0; i < int(mode.size()/3); i++){
        if(mode[i] == 0){
            bool ifremained = false;
            for (int k = 0; k < remain_idx.size(); k++){
                if (i == remain_idx[k]){
                    ifremained = true;
                }
            }
            if(!ifremained){
                envs.push_back(pts[i]);
                envs.back().d = 0.042;
            }
        }
    }
    for (int k = 0; k < remain_idx.size(); k++){
        envs.push_back(pts[remain_idx[k]]);
    }
    return VelocityCorrection(envs);
}

PlanningWorld::PlanningWorld(){
    this->cons = new ContactConstraints(2);
}

PlanningWorld::PlanningWorld(WorldTemplate* w){
    this->world = w;
    this->cons = new ContactConstraints(2);
}

Vector6d PlanningWorld::QPVelocity(const Vector6d& v_goal, const std::vector<ContactPoint>& mnps, 
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o, double wt, double wa){


    MatrixXd A_env;
    MatrixXd G_env;
    VectorXd b_env;
    VectorXd h_env;

    this->cons->ModeConstraints(envs, env_mode, this->mu_env, f_ext_o, &A_env, &b_env, &G_env, &h_env);
    // std::cout << "A_env:\n" << A_env << std::endl;
    // Vector6d vv;
    // vv << -1,1,0,0,0,1;
    // std::cout << A_env.block(0,0,6,6)*vv << std::endl;

    MatrixXd A_mnp;
    MatrixXd G_mnp;
    VectorXd b_mnp;
    VectorXd h_mnp;

    VectorXi mmode;
    mmode.resize(3*mnps.size());
    mmode.setZero(); // all fixed manipulator contacts

    this->cons->ModeConstraints(mnps, mmode, this->mu_mnp, f_ext_o, &A_mnp, &b_mnp, &G_mnp, &h_mnp);

    MatrixXd A;
    MatrixXd G;
    VectorXd b;
    VectorXd h;

    mergeManipulatorandEnvironmentConstraints_relax(A_mnp, b_mnp, G_mnp, h_mnp, A_env, b_env, G_env, h_env, &A, &b, &G, &h);

    int n_var = A.cols();

    MatrixXd P(n_var, n_var);
    VectorXd p(n_var);
    p.setZero();
    P.setIdentity();
    P.block(0,0,3,3) = wt*P.block(0,0,3,3);
    P.block(3,3,3,3) = wa*P.block(3,3,3,3);
    P.block(6,6,n_var-6, n_var-6) = 0.01*P.block(6,6,n_var-6, n_var-6);
    p.block(0,0,6,1) = -v_goal;

    VectorXd x(n_var);
    x.setZero();

    Vector6d x_v;

    // std::cout << "f: " << 
    double f = solve_quadprog(P, p, A.transpose(), -b,  G.transpose(), -h, x);
    // std::cout << "P \n" << P << std::endl; 
    // std::cout << "p \n" << p << std::endl; 
    // std::cout << "A.transpose() \n" << A.transpose() << std::endl; 
    // std::cout << "-b \n" << -b.transpose() << std::endl; 
    // std::cout << "-h \n" << -h.transpose() << std::endl; 
    // std::cout << "G.transpose() \n" << G.transpose() << std::endl; 
    // std::cout << "x: \n" << x << std::endl;
 
    
    // bool ifconstrained = ifConstraintsSatisfied(x, A, b, G,h);
    if (std::isinf(f)){ // if fail to solve the problem
        x_v.setZero();
    // } else if (!ifConstraintsSatisfied(x, A, b, G,h)) {
    //     std::cout << " Constraints not satisfied for qp! " << std::endl;
    //     x_v.setZero();
    } else {
        x_v = x.block(0,0,6,1);
    }
    
    return x_v;
    
}

Vector6d PlanningWorld::QPVelocity_Quasidynamic(const Vector6d& v_goal, const std::vector<ContactPoint>& mnps, 
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o, double wt, double wa,
    double h_time){


    MatrixXd A_env;
    MatrixXd G_env;
    VectorXd b_env;
    VectorXd h_env;

    this->cons->ModeConstraints(envs, env_mode, this->mu_env, f_ext_o, &A_env, &b_env, &G_env, &h_env);
    // std::cout << "A_env:\n" << A_env << std::endl;
    // Vector6d vv;
    // vv << -1,1,0,0,0,1;
    // std::cout << A_env.block(0,0,6,6)*vv << std::endl;

    MatrixXd A_mnp;
    MatrixXd G_mnp;
    VectorXd b_mnp;
    VectorXd h_mnp;

    VectorXi mmode;
    mmode.resize(3*mnps.size());
    mmode.setZero(); // all fixed manipulator contacts

    this->cons->ModeConstraints(mnps, mmode, this->mu_mnp, f_ext_o, &A_mnp, &b_mnp, &G_mnp, &h_mnp);

    MatrixXd A;
    MatrixXd G;
    VectorXd b;
    VectorXd h;

    mergeManipulatorandEnvironmentConstraints_relax(A_mnp, b_mnp, G_mnp, h_mnp, A_env, b_env, G_env, h_env, &A, &b, &G, &h);

    // add quasidynamic condition

    int n_var = A.cols();

    int G_rows = G.rows();
    if(G_rows < 12){

        G.conservativeResize(12 + G_rows, n_var);
        h.conservativeResize(12  + G_rows);

        G.block(G_rows, 0, 12, n_var).setZero();
        h.block(G_rows, 0, 12, 1).setZero();
        
        VectorXd sigma = VectorXd::Constant(f_ext_o.size(), 1e-4);
        h.block(G_rows, 0 , 6, 1) = -f_ext_o - sigma;
        h.block(G_rows + 6, 0 , 6, 1) = -(-f_ext_o + sigma);
    }
    
    G.block(G.rows()-12,0,6,6) = -this->object_inertia*1/(h_time); // MI: object inertia
    G.block(G.rows()-6,0,6,6) = this->object_inertia*1/(h_time);
    

    

    MatrixXd P(n_var, n_var);
    VectorXd p(n_var);
    p.setZero();
    P.setIdentity();
    P.block(0,0,3,3) = wt*P.block(0,0,3,3);
    P.block(3,3,3,3) = wa*P.block(3,3,3,3);
    P.block(6,6,n_var-6, n_var-6) = 0.01*P.block(6,6,n_var-6, n_var-6);
    p.block(0,0,6,1) = -v_goal;

    VectorXd x(n_var);
    x.setZero();

    Vector6d x_v;

    // std::cout << "f: " << 
    
    double f = solve_quadprog(P, p, A.transpose(), -b,  G.transpose(), -h, x);

    if (std::isinf(f)){ // if fail to solve the problem
        x_v.setZero();
    } else if (!ifConstraintsSatisfied(x, A, b, G, h)) {
        std::cout << " Constraints not satisfied for qp! " << std::endl;
        x_v.setZero();
        // x_v = x.block(0,0,6,1);
    } else {
        x_v = x.block(0,0,6,1);
    }
    return x_v;
    
}


bool PlanningWorld::ForwardIntegration(const Vector7d& x_start, const Vector7d& x_goal, 
    const VectorXd& mnp_config, const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_, const Vector6d& f_ext_w,
    double wt, double wa, std::vector<Vector7d>* path, bool ifquasidynamic){
    
    std::cout<<"Forward integration from: " << x_start.transpose() << std::endl;
    std::cout << "mnp_config: " << mnp_config.transpose() << std::endl;

    double thr = 1e-4;
    // double h = 0.02;
    // int max_counter = 150;
    double h = 0.04;
    // double h = 0.08;
    int max_counter = 150;

    bool ifmanipulatorcollide = false;

    Vector7d x = x_start;
    VectorXi env_mode = env_mode_;

    path->push_back(x);

    std::vector<ContactPoint> envs;
    std::vector<ContactPoint> mnps;
    std::vector<ContactPoint> fingertips;
    this->world->getRobot()->getFingertipsOnObject(mnp_config, x, &fingertips);
    this->world->getRobot()->Fingertips2PointContacts(fingertips, &mnps);
    // this->world->getObjectContacts(&envs, x_new)(x, &envs);
    envs = envs_;

    std::vector<ContactPoint> envs_pre;
    envs_pre = envs;

    Vector6d v_b_pre;
    v_b_pre.setZero();

    // if(envs.size()!=int(env_mode.size()/3)){
    //     std::cout << "Environment contacts and modes are not consistent" << std::endl;
    //     if(envs_.size()==2){
    //         std::vector<ContactPoint> envs_update;
    //         if(same_line_update(envs_, envs, &envs_update)){
    //             envs = envs_update;
    //         } else {
    //             std::cout << "failed to do same line update!" << std::endl;
    //             return true;
    //         }
    //     } else {
    //         deleteExtraContacts(envs_, envs);
    //         if (envs.size()!=envs_.size()){
    //             std::cout << "failed to delete extra contacts!" << std::endl;
    //             return true;
    //         }
    //     }
    // }
    int counter;
    int delete_c=0;

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
        if (ifquasidynamic){
            v_b = QPVelocity_Quasidynamic(v_star, mnps, envs, env_mode, f_ext_o, wt, wa, 1);
        } else {
            v_b = QPVelocity(v_star, mnps, envs, env_mode, f_ext_o, wt, wa);
        }

        // Vector6d v_quasidynamic = QPVelocity_Quasidynamic(v_star, mnps, envs, env_mode, f_ext_o, wt, wa, 1);
        // Vector6d v_quasistatic =  QPVelocity(v_star, mnps, envs, env_mode, f_ext_o, wt, wa);
        // double d2 = (v_quasidynamic - v_quasistatic).norm();
        
        // Vector6d v_b_2 = EnvironmentConstrainedVelocity(v_star, envs, env_mode);
        // if (!isForceBalance(mnps, envs, env_mode, f_ext_o)){
        //     v_b_2.setZero();
        // }
        
        if (v_b.norm() < thr ){
            std::cout << "v_b < thr : " << v_b.transpose() << std::endl;
            break;
        }

        if ((v_b_pre.transpose()*v_b)[0] < -1e-5){
            printf("v_b back and forth. \n");
            break;
        }
        
        steer_velocity(v_b, h, this->charac_len);

        // integrate v
        Vector7d x_new = SE32pose(T*se32SE3(v_b));

        // check manipulator collision (break may happen)
        this->world->updateObjectPose(x_new);

        if(this->world->getRobot()->ifIKsolution(mnp_config, x_new)){
            ifmanipulatorcollide = this->world->isRobotCollide(mnp_config);
        } else {
            ifmanipulatorcollide = true;
        }

        // std::cout << "Manipulator collision:" << ifmanipulatorcollide << std::endl;
        if (ifmanipulatorcollide){
                printf("Manipulator collide!\n");
                break;
        }

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

        if(counter == max_counter - 1){
            printf("Reach the end.\n");
        }
    }

    std::cout<<"counter:" << counter << " x: " << x.transpose() << std::endl;

    return ifmanipulatorcollide;
}


bool PlanningWorld::ForwardIntegration_ObjectPath(const Vector7d& x_start, const Vector7d& x_goal, 
    const std::vector<ContactPoint>& envs_, const VectorXi& env_mode_,
    double wt, double wa, std::vector<Vector7d>* path){
    
    std::cout<<"Forward integration from: " << x_start.transpose() << std::endl;

    double thr = 1e-4;

    double h = 0.04;

    int max_counter = 150;

    Vector7d x = x_start;
    VectorXi env_mode = env_mode_;

    path->push_back(x);

    std::vector<ContactPoint> envs;
    envs = envs_;

    std::vector<ContactPoint> envs_pre;
    envs_pre = envs;

    Vector6d v_b_pre;
    v_b_pre.setZero();

    int counter;
    int delete_c=0;

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

        Vector6d v_b = EnvironmentConstrainedVelocity(v_star, envs, env_mode);
        
        if (v_b.norm() < thr ){
            std::cout << "v_b < thr : " << v_b.transpose() << std::endl;
            break;
        }

        if ((v_b_pre.transpose()*v_b)[0] < -1e-5){
            printf("v_b back and forth. \n");
            break;
        }
        
        steer_velocity(v_b, h, this->charac_len);

        // integrate v
        Vector7d x_new = SE32pose(T*se32SE3(v_b));

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

        if(counter == max_counter - 1){
            printf("Reach the end.\n");
        }
    }

    std::cout<<"counter:" << counter << " x: " << x.transpose() << std::endl;

    return true;
}


Vector6d PlanningWorld::EnvironmentConstrainedVelocity(const Vector6d& v_goal, 
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode){
    
    MatrixXd A_env;
    MatrixXd G_env;
    VectorXd b_env;
    VectorXd h_env;

    this->cons->ModeConstraints(envs, env_mode, this->mu_env, Vector6d::Zero(), &A_env, &b_env, &G_env, &h_env);

    MatrixXd A;
    MatrixXd G;
    VectorXd b;
    VectorXd h;

    // std::cout << "A_env\n" << A_env.block(0,0,A_env.rows(),6) << std::endl;
    // std::cout << "G_env\n" << G_env.block(0,0,G_env.rows(),6) << std::endl;

    deleteZeroRows(A_env.block(0,0,A_env.rows(),6), b_env, &A, &b);
    deleteZeroRows(G_env.block(0,0,G_env.rows(),6), h_env, &G, &h);

    // std::cout << "A\n" << A << std::endl;
    // std::cout << "G\n" << G << std::endl;

    mergeDependentRows(A, b, &A, &b);

    int n_var = 6;

    MatrixXd P(6, 6);
    VectorXd p(6);
    P.setIdentity();
    p = -v_goal;

    VectorXd x(n_var);
    //x.setZero();
    x = v_goal;

    if(A.rows() > n_var){
        FullPivLU<MatrixXd> lu_decomp(A.transpose());

        if(lu_decomp.rank() >= n_var){
        // if A fully constrainted the velocity
            x.setZero();
            // double f = solve_quadprog(P, p, A.transpose(), -b,  G.transpose(), -h, x);
            return x;
        } else {
            A = (lu_decomp.image(A.transpose())).transpose();
            b = VectorXd::Zero(A.rows());
        }
    }

    // std::cout << "f: " << 
    // std::cout << "G \n" << P << std::endl; 
    // std::cout << "g0 \n" << p << std::endl; 
    // std::cout << "CE \n" << A.transpose() << std::endl; 
    // std::cout << "ce \n" << -b.transpose() << std::endl; 
    // std::cout << "CI \n" << G.transpose() << std::endl; 
    // std::cout << "ci \n" << -h.transpose() << std::endl; 
    double f = solve_quadprog(P, p, A.transpose(), -b,  G.transpose(), -h, x);
    // std::cout << "x: \n" << x << std::endl;

    // if (std::isinf(f) || (!ifConstraintsSatisfied(x, A, b, G, h)) || f > 0.0){ // if fail to solve the problem
    //     x.setZero();
    // } 
    if (std::isinf(f) || (!ifConstraintsSatisfied(x, A, b, G, h))){ // if fail to solve the problem
        x.setZero();
    } 

    // if(f > 0.0){
    //     std::cout << "solve_quadprog error" << std::endl;
    // }
    
    return x;

}

bool PlanningWorld::isForceBalance(const std::vector<ContactPoint>& mnps, 
    const std::vector<ContactPoint>& envs, const VectorXi& env_mode, const Vector6d& f_ext_o){
    
    if (mnps.size() + envs.size() == 0){
        return false;
    }
    MatrixXd A_env;
    MatrixXd G_env;
    VectorXd b_env;
    VectorXd h_env;

    this->cons->ModeConstraints(envs, env_mode, this->mu_env, f_ext_o, &A_env, &b_env, &G_env, &h_env);
    // std::cout << "A_env:\n" << A_env << std::endl;
    // Vector6d vv;
    // vv << -1,1,0,0,0,1;
    // std::cout << A_env.block(0,0,6,6)*vv << std::endl;

    MatrixXd A_mnp;
    MatrixXd G_mnp;
    VectorXd b_mnp;
    VectorXd h_mnp;

    VectorXi mmode;
    mmode.resize(3*mnps.size());
    mmode.setZero(); // all fixed manipulator contacts

    this->cons->ModeConstraints(mnps, mmode, this->mu_mnp, f_ext_o, &A_mnp, &b_mnp, &G_mnp, &h_mnp);

    // std::cout <<"A_mnp" << A_mnp << std::endl;

    MatrixXd A;
    MatrixXd G;
    VectorXd b;
    VectorXd h;

    mergeManipulatorandEnvironmentConstraints(A_mnp, b_mnp, G_mnp, h_mnp, A_env, b_env, G_env, h_env, &A, &b, &G, &h);
    
    int n_var = A.cols() - 6;
    if(n_var == 0){
        return false;
    }
    A = A.block(0, 6, A.rows(), n_var);
    G = G.block(0, 6, G.rows(), n_var);

    deleteZeroRows(A, b, &A, &b);
    deleteZeroRows(G, h, &G, &h);

    VectorXd xl;// = VectorXd::Constant(n_var, std::numeric_limits<double>::min());
    VectorXd xu;// = VectorXd::Constant(n_var, std::numeric_limits<double>::max());
    
    VectorXd C = VectorXd::Constant(n_var, 0);
    VectorXd x(n_var); 
    x.setZero();
    double optimal_cost;

    // bool result = lp(C, -G, -h, A, b, xl, xu, &x, &optimal_cost);

    VectorXd sigma = VectorXd::Constant(b.size(), 1e-4);
    MatrixXd A_relax;
    MatrixXd G_relax(G.rows()+2*A.rows(), n_var);
    VectorXd b_relax;
    VectorXd h_relax(G.rows()+2*A.rows());
    h_relax.setZero();
    G_relax.block(0,0,G.rows(),n_var) = G;
    G_relax.block(G.rows(),0,A.rows(),n_var) = A;
    G_relax.block(G.rows()+A.rows(),0,A.rows(),n_var) = -A;
    h_relax.block(G.rows(),0,A.rows(),1) = b - sigma;
    h_relax.block(G.rows()+A.rows(),0,A.rows(),1) = -(b + sigma);

    // std::cout << "G_relax\n" << G_relax << std::endl;
    // std::cout << "h_relax\n" << h_relax << std::endl;

    bool result_relax = lp(C, -G_relax, -h_relax, A_relax, b_relax, xl, xu, &x, &optimal_cost);
    // if(result){
    //     std::cout << "A\n" << A << std::endl;
    //     std::cout << "b\n" << b << std::endl;
    //     std::cout << "G\n" << G << std::endl;
    //     std::cout << "h\n" << h << std::endl;
    // }

    return result_relax;
}