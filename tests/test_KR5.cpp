// #include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartKR5.h"

#include "cmgmp/worlds/DartWorld.h"
// #include "cmgmp/worlds/PlanningWorld.h"

#include <ctime>
#include <iostream>
#include <fstream> 

int main(int argc, char* argv[])
{
    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", 0.1*Vector3d(1,1,1));

    SkeletonPtr env1 = createFixedBox("ground", 0.1*2*Vector3d(10,10,1), 0.1*Vector3d(0,0,-10));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);

    Vector7d x;
    x << 0,0,0,0,0,0,1;
    my_world->updateObjectPose(x);
    

    DartKR5* rpt = new DartKR5();
    my_world->addRobot(rpt);

    // Position its base in a reasonable way
    Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
    tf.translation() = Eigen::Vector3d(0, 0.0, 0.0);
    rpt->bodies[0]->getJoint(0)->setTransformFromParentBodyNode(tf);

    // std::vector<Vector3d> pts;
    // rpt->computeWorkspace(&pts);

    // std::vector<ContactPoint> pts_world;
    VectorXd joint_v;
    VectorXd p(6);
    // p << 0.2,0.2,0.3,0,0,-1;
    // p << 0.7,0,0.025,0,0,-1;
    // p << 0.47,0,0.023*1.5,0,0,-1;
    // p << 0.482697,-0.025495,0.0082869,0,0,-1;
    p << 0.4,0,0.4,0,0,-1;





    bool iksolved = rpt->inverseKinematics(p, x, joint_v);
    std::cout << "if ik solved " << iksolved  << ", joint values: " << joint_v.transpose() << std::endl;
    std::cout << "if robot collide" << my_world->isRobotCollide(p) << std::endl;
    // std::cout << "if robot collide: " << my_world->isRobotCollide(r); 
    my_world->startWindow(&argc, argv);
}