#include "cmgmp/search/rrtplanner.h"
#include "cmgmp/utilities/sample_grasp.h"
#include "cmgmp/contacts/contact_kinematics.h"

#include "cmgmp/contacts/contact_utils.h"

#include "cmgmp/manipulators/DartDDHand.h"

#include "cmgmp/worlds/DartWorld.h"

#include <ctime>
#include <iostream>
#include <fstream> 

int main(int argc, char* argv[])
{
    DartWorld* my_world = new DartWorld();

    
    SkeletonPtr object = createFreeBox("box_object", 0.1*Vector3d(1,1,1));

    SkeletonPtr env1 = createFixedBox("ground", 0.1*2*Vector3d(10,10,1), 0.1*Vector3d(0,0,-1));

    my_world->addObject(object);
    my_world->addEnvironmentComponent(env1);

    Vector7d x;
    x << 0,0,0.05,0,0,0,1;
    my_world->updateObjectPose(x);
    

    DartDDHand* rpt = new DartDDHand(L_FINGER);
    my_world->addRobot(rpt);

    // Position its base in a reasonable way
    Vector7d handx;
    handx << -0.2,0,55E-3 + 28E-3,0,0,0,1;
    rpt->setHandFrameTransform(handx);


    for (double xx = 0.0; xx < 55E-3; xx += 1E-3){
        for (double zz = 0.0; zz < 90E-3; zz += 1E-3){
            bool if_hand_ik = rpt->hand_ik->is_inverse_kinematics(Vector2d(-xx, -zz), Vector2d(xx, -zz));
            if (if_hand_ik){

                std::cout << "hand ik true: [" << xx << ", " << zz << "]," << std::endl;
            }
        }
    }

    // bool if_hand_ik = rpt->hand_ik->is_inverse_kinematics(Vector2d(xx, zz), Vector2d(-xx, zz));
    // std::cout << "if hand ik " << if_hand_ik << std::endl;

    my_world->startWindow(&argc, argv);
}