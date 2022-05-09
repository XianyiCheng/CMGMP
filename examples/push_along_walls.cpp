#include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartPointManipulator.h"

#include "cmgmp/worlds/DartWorld.h"
// #include "cmgmp/worlds/PlanningWorld.h"

#include <ctime>
#include <iostream>
#include <fstream>  

#define VISUALIZE_SG  0
#define VISUALIZE_PTS  1
#define DO_SEARCH  2

typedef Matrix<double, 12, 1> Vector12d;

void push_along_walls(Vector7d& x_start, Vector7d& x_goal, double& goal_thr, 
    double& wa, double& wt, double &eps_trans, double& eps_angle,
    Vector3d& x_ub, Vector3d& x_lb, PlanningWorld* pw, std::vector<ContactPoint>* surface){
    
    // parameters 
    pw->mu_env = 0.4;
    pw->mu_mnp = 0.8;

    x_ub << 5, 11, 0;
    x_lb << -5, 1, 0;
    
    wa = 0.3;
    wt = 1;

    eps_trans = 10;
    eps_angle = 3.14*90/180;
    goal_thr = 3.14*5/180;

    // load object surface discretization;
    surfacepoints_from_file("../data/grasp_sampling/box_halflength_1.csv", surface);
    for (int i = 0; i < surface->size(); i++){
        surface->at(i).p = 0.1*surface->at(i).p;
    }

    // x_start<< 0.1*5, 0, 0.1*1, 0, 0, 0, 1;
    x_start<< 0.5, 0, 0.1, 0., 0, 0, 1;
    x_goal<< 0.1*-3, 0, 1.1, 0, 0, 0, 1;
    // x_start<< 5, 1, 0, 0, 1, 0, 0;
    // x_goal<< 0.1*-3, 0, 0.1*11, 0, 0, 0, 1;
    // x_goal<< -3, 11, 0, 0, 1, 0, 0;

    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", 0.1*Vector3d(2,2,2));

    SkeletonPtr env1 = createFixedBox("ground", 0.1*2*Vector3d(10,10,5), 0.1*Vector3d(10,0,-5));
    SkeletonPtr env2 = createFixedBox("wall", 0.1*2*Vector3d(5,10,10), 0.1*Vector3d(-5,0,0));

    // SkeletonPtr env1 = createFixedBox("ground", 2*Vector3d(10,5,10), Vector3d(10,-5,0));
    // SkeletonPtr env2 = createFixedBox("wall", 2*Vector3d(5,10,10), Vector3d(-5,0,0));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);
    my_world->addEnvironmentComponent(env2);
    

    DartPointManipulator* rpt = new DartPointManipulator(2, 0.1);
    my_world->addRobot(rpt);

    pw->world = my_world;
    std::cout <<pw->world->getRobot()->getNumberofFingertips() << std::endl;
  
}

int main(int argc, char* argv[]){

    RRTPlannerOptions options;
    options.goal_biased_prob = 1;
    options.max_samples = 10;
    options.sampleSO3 = false;
    options.rotation_sample_axis << 0,0,1;

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
    f_w << 0,0,-1,0,0,0;

    push_along_walls(x_start, x_goal, goal_thr, wa, wt, eps_trans, eps_angle, x_ub, x_lb, &pw, &surface);

    
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