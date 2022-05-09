#include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartDDHandFreeY.h"

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
double hand_z = 55E-3 + 28E-3 + 0.5*34E-3;

void save_hand_frame_traj(std::string filename, RRTPlanner & planner, const std::vector<int>& path){

    DartDDHandFreeY* hand = new DartDDHandFreeY(L_FINGER);
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


    MatrixXd mm(traj.size(), traj[0].size());

    for (int i = 0; i < traj.size(); i++){
        mm.row(i) = traj[i].transpose();
    }
    saveData(filename, mm);
}


void bookshelf(Vector7d& x_start, Vector7d& x_goal, double& goal_thr, 
    double& wa, double& wt, double &eps_trans, double& eps_angle,
    Vector3d& x_ub, Vector3d& x_lb, PlanningWorld* pw, std::vector<ContactPoint>* surface){
    
    // parameters 
    pw->mu_env = 0.8;
    pw->mu_mnp = 0.8;
    pw->charac_len = 50.0;

    x_ub << 0, 0, 17;
    x_lb << 0, -34E-3, 34E-3;
    
    wa = 0.5;
    wt = 1;

    double scale = 0.01;


    eps_trans = scale*1;
    eps_angle = 3.14*15/180;
    goal_thr = 3.14*5/180;


    // load object surface discretization;
    surfacepoints_from_file("../data/grasp_sampling/box_lx_1_ly_4_lz_4.csv", surface);
    for (std::vector<ContactPoint>::iterator it = surface->begin(); it != surface->end(); ) 
    {    
        if(((it->p[2] > 1.7) && (abs(it->n[0] > 0.5))) || (it->p[1] > -0.5) || (it->p[1] < -1.99 && it->n[1] > 0.5)) 
            it = surface->erase(it);
            
        else{ 
            it->p = Vector3d(10E-3/1 * it->p[0], 34E-3/4 * it->p[1], 34E-3/4 * it->p[2]);
            ++it;
        }
    }

    x_start<< 0.0, 0.0 ,0.99*17E-3, 0, 0, 0, 1;
    // x_goal<< 0.0, -17E-3, 34E-3 , 0.25881904510252074, 0, 0, 0.965925826289;
    x_goal<< 0.0, -17E-3, 34E-3 ,0.17364817766693033, 0, 0, 0.984807753012208;
    // x_goal<< 0.0, -34E-3, 34E-3 , 0.7071, 0, 0, 0.7071;
    // x_goal<< 0.0, -34E-3, 34E-3 , 0.5, 0, 0, 0.8660;



    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", Vector3d(10E-3,34E-3,34E-3));

    SkeletonPtr env0 = createFixedBox("ground", Vector3d(0.5, 0.5, 0.2), Vector3d(0,0,-0.1));
    SkeletonPtr env1 = createFixedBox("book1", Vector3d(10E-3,34E-3,34E-3), Vector3d(10E-3 + 1E-3,0,17E-3));
    SkeletonPtr env2 = createFixedBox("book2", Vector3d(10E-3,34E-3,34E-3), Vector3d(-(10E-3 + 1E-3),0,17E-3));
    SkeletonPtr env3 = createFixedBox("back", Vector3d(3*10E-3 + 3E-3, 10E-3, 34E-3), Vector3d(0,17E-3 + 5E-3 + 2E-3, 17E-3));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env0);
    my_world->addEnvironmentComponent(env1);
    my_world->addEnvironmentComponent(env2);
    my_world->addEnvironmentComponent(env3);
    

    DartDDHandFreeY* rpt = new DartDDHandFreeY(I_FINGER);
    Vector7d handx;
    handx << hand_x,hand_y,hand_z,0,0,0,1;
    
    rpt->setHandFrameTransform(handx);

    my_world->addRobot(rpt);

    pw->world = my_world;
  
}

int main(int argc, char* argv[]){

    RRTPlannerOptions options;
    options.goal_biased_prob = 0.6;
    options.max_samples = 200;
    options.sampleSO3 = true;

    // options.sampleSO3 = false;
    // options.rotation_sample_axis << 1,0,0;

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

    bookshelf(x_start, x_goal, goal_thr, wa, wt, eps_trans, eps_angle, x_ub, x_lb, &pw, &surface);

    
    if (test_option == VISUALIZE_SG){
        std::vector<Vector7d> pp;
        pp.push_back(x_start);
        // pp.push_back(x_goal);

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

        planner.Search(options, x_goal, goal_thr, &path, success, t);

        planner.VisualizePath(path);
        
        save_hand_frame_traj("ddhand_bookshelf.csv", planner, path);

    }

    pw.world->startWindow(&argc, argv);

}
