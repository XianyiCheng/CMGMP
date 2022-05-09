#include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartDDHand.h"

#include "cmgmp/worlds/DartWorld.h"
// #include "cmgmp/worlds/PlanningWorld.h"

#include <ctime>
#include <iostream>
#include <fstream>  

#define VISUALIZE_SG  0
#define VISUALIZE_PTS  1
#define DO_SEARCH  2

double hand_x = 0;
double hand_y = 0.0;
double hand_z = 55E-3 + 28E-3;

void save_hand_frame_traj(std::string filename, RRTPlanner & planner, const std::vector<int>& path){

    DartDDHand* hand = new DartDDHand(L_FINGER);
    Vector7d handx;
    handx << hand_x,hand_y,hand_z,0,0,0,1;

    hand->setHandFrameTransform(handx);

    std::vector<VectorXd> traj;

    for(int i = 0; i < path.size(); i++){
        int k = path[i];

        if(k==0){continue;}
        std::vector<Vector7d> path = planner.T->edges[planner.T->nodes[k].edge].path;
        VectorXd mnp_config = planner.T->nodes[k].manipulator_config;
        for (auto x: path){
            VectorXd rc = hand->getHandFrameConfig(mnp_config, x);
            traj.push_back(rc);
        }
    }

    MatrixXd mm(traj.size(), 4);

    for (int i = 0; i < traj.size(); i++){
        mm.row(i) = traj[i].transpose();
    }
    saveData(filename, mm);
}


void flip_n_pinch(Vector7d& x_start, Vector7d& x_goal, double& goal_thr, 
    double& wa, double& wt, double &eps_trans, double& eps_angle,
    Vector3d& x_ub, Vector3d& x_lb, PlanningWorld* pw, std::vector<ContactPoint>* surface){
    
    // parameters 
    pw->mu_env = 0.4;
    pw->mu_mnp = 0.8;
    pw->charac_len = 10.0;

    x_ub << 0.2, 0.02, 0;
    x_lb << -0.2, 0.005, 0;
    
    wa = 1;
    wt = 0.05;

    double scale = 0.01;


    eps_trans = scale*1;
    eps_angle = 3.14*120/180;
    goal_thr = 3.14*30/180;


    // load object surface discretization;
    surfacepoints_from_file("../data/grasp_sampling/box_lx_1_ly_4_lz_4.csv", surface);
    for (std::vector<ContactPoint>::iterator it = surface->begin(); it != surface->end(); ) 
    {
        if((abs(it->p[1]) > 0.25)) 
            it = surface->erase(it);
        else {
            // it->p = scale*it->p;
            it->p = Vector3d(10E-3/1 * it->p[0], 34E-3/4 * it->p[1], 34E-3/4 * it->p[2]);
            ++it;
        }

    }

    x_start<< 0.0, 0,0.4999*scale, 0, 0.7071, 0, 0.7071;
    
    // x_goal<< 0, 0,0.4999*scale, 0, -0.7071, 0, 0.7071;
    // x_goal<< 0.0, 0,(1 + 1.7)*0.5*scale , 0, -0.5, 0, 0.8660;
    x_goal<< 0.0, 0,(1 + 1.7)*0.5*scale , 0, -0.5, 0, 0.8660;

    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", scale*Vector3d(1,4,4));

    SkeletonPtr env1 = createFixedBox("ground", Vector3d(0.5, 0.5, 0.2), Vector3d(0,0,-0.1));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);
    

    DartDDHand* rpt = new DartDDHand(L_FINGER);
    Vector7d handx;
    handx << hand_x,hand_y,hand_z,0,0,0,1;

    rpt->setBoxConstraints(FINGER_1, Vector3d(0.1, 30E-3, -53E-3), Vector3d(-0.1, -30E-3, -120E-3));
    rpt->setBoxConstraints(FINGER_0, Vector3d(0.1, 30E-3, -53E-3), Vector3d(-0.1, -30E-3, -120E-3));
    
    rpt->setHandFrameTransform(handx);

    my_world->addRobot(rpt);

    pw->world = my_world;
    // std::cout <<pw->world->getRobot()->getNumberofFingertips() << std::endl;

    Matrix6d oi;
    oi << 1,0,0,0,0,0,
          0,1,0,0,0,0,
          0,0,1,0,0,0,
          0,0,0,1.0/6,0,0,
          0,0,0,0,1.0/6,0,
          0,0,0,0,0,1.0/6;
    pw->object_inertia = 0.05*oi; 
  
}

int main(int argc, char* argv[]){

    RRTPlannerOptions options;
    options.goal_biased_prob = 1;
    options.max_samples = 50;
    options.sampleSO3 = true;
    // options.rotation_sample_axis << 0,0,1;

    int test_option = DO_SEARCH;

    Vector7d x_start;
    Vector7d x_goal;
    PlanningWorld pw;
    double wa;
    double wt;
    Vector3d x_ub;
    Vector3d x_lb;
    double eps_trans;
    double eps_angle;
    double goal_thr;
    std::vector<ContactPoint> surface;
    Vector6d f_w;
    f_w << 0,0,-0.05,0,0,0;

    flip_n_pinch(x_start, x_goal, goal_thr, wa, wt, eps_trans, eps_angle, x_ub, x_lb, &pw, &surface);

    
    if (test_option == VISUALIZE_SG){
        std::vector<Vector7d> pp;
        pp.push_back(x_start);
        pp.push_back(x_goal);

        pw.world->setObjectTrajectory(pp);

    } 
    else if (test_option == VISUALIZE_PTS){
        Vector7d x;
        x << 0,0,0,0,0,0,1;
        pw.world->setSurfacePoints(surface); // TODO: draw contact points and their normals
    } 
    else if( test_option == DO_SEARCH)
    {

        RRTPlanner planner(&pw);

        planner.Initialize(x_lb, x_ub, x_start, f_w, surface, wa, wt, eps_trans, eps_angle);
        
        // planner.SetInitialManpulatorLocations(mnps_config);

        std::vector<int> path;
        
        double t;
        bool success;

        planner.Search(options, x_goal, goal_thr, &path, success, t, true);

        planner.VisualizePath(path);
        
        save_hand_frame_traj("ddhand_flip_n_pinch.csv", planner, path);

    }

    pw.world->startWindow(&argc, argv);

}
