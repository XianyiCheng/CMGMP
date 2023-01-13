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

#define VISUALIZE_SG 0
#define VISUALIZE_PTS 1
#define DO_SEARCH 2

typedef Matrix<double, 12, 1> Vector12d;

void card(Vector7d &x_start, Vector7d &x_goal, double &goal_thr,
          double &wa, double &wt, double &eps_trans, double &eps_angle,
          Vector3d &x_ub, Vector3d &x_lb, PlanningWorld *pw, std::vector<ContactPoint> *surface)
{
    double box_length = 2.0;
    double box_height = 0.5;

    DartWorld *my_world = new DartWorld();

    SkeletonPtr object =
        createFreeBox("box_object", Vector3d(box_length, box_length, box_height));
    SkeletonPtr env1 =
        createFixedBox("ground", Vector3d(4, 4, 0.2), Vector3d(0, 0, -0.1));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);

    int n_robot_contacts = 2;
    DartPointManipulator *rpt =
        new DartPointManipulator(n_robot_contacts, 0.15);
    my_world->addRobot(rpt);
    rpt->is_patch_contact = true;

    pw->world = my_world;

    // parameters
    pw->mu_env = 0.2;
    pw->mu_mnp = 0.7;

    x_start << 0, 0, box_height / 2 * 0.9999, 0, 0, 0, 1;

    // goal: rotate around y axis for 90 degrees
    x_goal << 2.5, 0, box_length / 2 * 1.5, 0, -0.7071, 0, 0.7071;

    x_ub << 5, 1, box_height;
    x_lb << -1, -1, 0.0;

    wa = 0.4;
    wt = 1.0;

    eps_trans = 0.75;
    eps_angle = 3.14 * 35 / 180;
    goal_thr = box_length * 3.14 * 10 / 180;

    // load object surface discretization;
    // MatrixXd data = load_points_from_csv(std::string(SRC_DIR) + "/data/grasp_sampling/peg_out_of_hole_surface_points.csv");
    // int N = data.rows();
    // for (int i = 0; i < N; ++i)  {
    //     ContactPoint pt((data.block(i,0,1,3)).transpose(), data.block(i,3,1,3).transpose(), 0);
    //     surface->push_back(pt);
    // }
    surfacepoints_from_file(std::string(SRC_DIR) + +"/data/grasp_sampling/box_halflength_1.csv", surface);
    for (std::vector<ContactPoint>::iterator it = surface->begin(); it != surface->end();)
    {
        if (abs(it->n[2]) < 0.8)
        {
            it = surface->erase(it);
        }
        else
        {
            ++it;
        }
    }
    for (int i = 0; i < surface->size(); i++)
    {
        surface->at(i).p = Vector3d(surface->at(i).p(0) * box_length / 2, surface->at(i).p(1) * box_length / 2, surface->at(i).p(2) * box_height / 2);
    }
}

int main(int argc, char *argv[])
{

    RRTPlannerOptions options;
    options.goal_biased_prob = 0.7;
    options.max_samples = 500;
    // options.sampleSO3 = false;
    // options.rotation_sample_axis << 0, 0, 1;

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
    f_w << 0, 0, -0.1, 0, 0, 0;

    card(x_start, x_goal, goal_thr, wa, wt, eps_trans, eps_angle, x_ub, x_lb, &pw, &surface);

    if (test_option == VISUALIZE_SG)
    {
        std::vector<Vector7d> pp;
        pp.push_back(x_start);
        pp.push_back(x_goal);

        pw.world->setObjectTrajectory(pp);
    }
    else if (test_option == VISUALIZE_PTS)
    {
        Vector7d x;
        x << 0, 0, 0, 0, 0, 0, 1;
        pw.world->setSurfacePoints(surface); // TODO: draw contact points and their normals
    }
    else if (test_option == DO_SEARCH)
    {

        RRTPlanner planner(&pw);

        planner.Initialize(x_lb, x_ub, x_start, f_w, surface, wa, wt, eps_trans, eps_angle);

        // planner.SetInitialManpulatorLocations(mnps_config);

        std::vector<int> path;

        double t;
        bool success;

        planner.Search(options, x_goal, goal_thr, &path, success, t, false);
        double grasp_measure_charac_length = 1.0;
        planner.printResults(path, success, t, grasp_measure_charac_length);
        planner.VisualizePath(path);
        //
    }

    pw.world->startWindow(&argc, argv);
}