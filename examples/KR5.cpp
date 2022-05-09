#include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartKR5.h"

#include "cmgmp/worlds/DartWorld.h"
// #include "cmgmp/worlds/PlanningWorld.h"

#include <ctime>
#include <iostream>
#include <fstream> 

#define VISUALIZE_SG  0
#define VISUALIZE_PTS  1
#define DO_SEARCH  2

void pushing(Vector7d& x_start, Vector7d& x_goal, double& goal_thr, 
    double& wa, double& wt, double &eps_trans, double& eps_angle,
    Vector3d& x_ub, Vector3d& x_lb, PlanningWorld* pw, std::vector<ContactPoint>* surface){
    
    // parameters 
    pw->mu_env = 0.4;
    pw->mu_mnp = 0.8;

    x_ub << 1, 1, 0.1;
    x_lb << -1, -1, 0.05;
    
    wa = 0.3;
    wt = 10;

    eps_trans = 0.5;
    eps_angle = 3.14*90/180;
    goal_thr = 0.05*3.14*30/180;

    // load object surface discretization;
    surfacepoints_from_file("../data/grasp_sampling/box_halflength_1.csv", surface);
    for (int i = 0; i < surface->size(); i++){
        surface->at(i).p = 0.025*surface->at(i).p;
    }

    // x_start<< 0.1*5, 0, 0.1*1, 0, 0, 0, 1;
    x_start<< 0.48, 0, 0.025*0.9999, 0, 0, 0, 1;
    // x_start<< 5, 1, 0, 0, 1, 0, 0;
    x_goal<< 0.48, 0, 0.025*0.9999, 0, 0, 0, 1;
    // x_goal<< -3, 11, 0, 0, 1, 0, 0;

    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", 0.05*Vector3d(1,1,1));
    SkeletonPtr env1 = createFixedBox("ground", Vector3d(2, 2, 0.2), Vector3d(0,0,-0.1));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);
    
    DartKR5* rpt = new DartKR5();
    my_world->addRobot(rpt);

    pw->world = my_world;
    std::cout <<pw->world->getRobot()->getNumberofFingertips() << std::endl;

    Matrix6d oi;
    oi << 1,0,0,0,0,0,
          0,1,0,0,0,0,
          0,0,1,0,0,0,
          0,0,0,1.0/6,0,0,
          0,0,0,0,1.0/6,0,
          0,0,0,0,0,1.0/6;
    pw->object_inertia = oi; 
  
}

void pivoting(Vector7d& x_start, Vector7d& x_goal, double& goal_thr, 
    double& wa, double& wt, double &eps_trans, double& eps_angle,
    Vector3d& x_ub, Vector3d& x_lb, PlanningWorld* pw, std::vector<ContactPoint>* surface){
    
    // parameters 
    pw->mu_env = 0.4;
    pw->mu_mnp = 0.8;
    pw->charac_len = 3.14/0.025;

    x_ub << 0.6, 0.02, 0.03;
    x_lb << 0.45, -0.02, 0.02;
    
    wa = 1;
    wt = 0.1;

    eps_trans = 0.01;
    eps_angle = 3.14*120/180;
    goal_thr = 3.14*10/180;

    // load object surface discretization;
    surfacepoints_from_file("../data/grasp_sampling/box_halflength_1.csv", surface);
    // add more candidates of surface points
    for (int i = 0; i < surface->size(); i++){
        surface->at(i).p = 0.02*surface->at(i).p;
    }

    // x_start<< 0.1*5, 0, 0.1*1, 0, 0, 0, 1;
    x_start<< 0.5, 0,0.01999 + 0.1, 0, 0, 0, 1;
    // x_goal << 0.5, 0,0.01999 + 0.1,-0.500000,-0.500000,0.500000,0.500000;
    x_goal << 0.5, 0,0.01999 + 0.1,-0.500000,0.500000,-0.500000,0.500000;
    // x_start << 0.516811,  1.34679e-16,    0.0342939, -9.29332e-16,     0.272343,  2.21832e-16,       0.9622;
    // x_start<< 0.518809,  0.00153148,   0.0333371,  -0.0315828,    0.420204, -0.00161959,    0.906878;
    // x_start<< 5, 1, 0, 0, 1, 0, 0;
    // x_goal<< 0.5, 0, 0.01999 + 0.1, 0, 0.7071, 0, 0.7071;
    // x_goal<< -3, 11, 0, 0, 1, 0, 0;

    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", 0.04*Vector3d(1,1,1));
    SkeletonPtr env1 = createFixedBox("ground", Vector3d(2, 2, 0.2), Vector3d(0,0,0));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);
    
    DartKR5* rpt = new DartKR5();
    my_world->addRobot(rpt);

    pw->world = my_world;
    std::cout <<pw->world->getRobot()->getNumberofFingertips() << std::endl;

    Matrix6d oi;
    oi << 1,0,0,0,0,0,
          0,1,0,0,0,0,
          0,0,1,0,0,0,
          0,0,0,1.0/6,0,0,
          0,0,0,0,1.0/6,0,
          0,0,0,0,0,1.0/6;
    pw->object_inertia = 0.01*oi; 
  
}


int main(int argc, char* argv[])
{
   RRTPlannerOptions options;
    options.goal_biased_prob = 0.5;
    options.max_samples = 500;
    options.sampleSO3 = true;
    // options.rotation_sample_axis << 0,0,1;

    int test_option = VISUALIZE_SG;

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
    f_w << 0,0,-0.01,0,0,0;

    pivoting(x_start, x_goal, goal_thr, wa, wt, eps_trans, eps_angle, x_ub, x_lb, &pw, &surface);

    
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

    }

    pw.world->startWindow(&argc, argv);
}