#include "DartWorld.h"
#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "../contacts/contact_kinematics.h"
#endif 

#ifndef DART_UTILS
#define DART_UTILS
    #include "../dart_utils/dart_utils.h"
#endif

#include <dart/collision/bullet/BulletCollisionDetector.hpp>

void DartWorldWindow::timeStepping() 
{

    if (play_mode == PLAY_BACK){
        if (object_positions.size() > 0 && (robot_configs.size() == object_positions.size())){
            int N = 50;
            if (int(frameCount/N) >= object_positions.size()){
            frameCount = 0;
            }
            object->setPositions(object_positions[int(frameCount/N)]);
            robot->setConfig(robot_configs[int(frameCount/N)], pose6d_to_pose7d(object_positions[int(frameCount/N)]));
            
            frameCount ++;
            // SimWindow::timeStepping();
        }
    }
    if (play_mode == SHOW_OBJECT_ALL){
        if (frameCount >= object_positions.size()){
            frameCount = 0;
        }
        object->setPositions(object_positions[frameCount]);
        frameCount++;
    }
}

DartWorld::DartWorld()
{
    world = std::make_shared<World>();
    window = std::make_shared<DartWorldWindow>(world);

    // bullet_collision_detector = BulletCollisionDetector::create();
    // auto fcl_collision = FCLCollisionDetector::create();
    // auto b_collision = DARTCollisionDetector::create();
    
    // world->getConstraintSolver()->setCollisionDetector(DARTCollisionDetector::create());

    world->getConstraintSolver()->setCollisionDetector(BulletCollisionDetector::create());
}

void DartWorld::addEnvironmentComponent(const SkeletonPtr& env){
    world->addSkeleton(env);
    environment.push_back(env);

    if (environmentCollisionGroup == 0){
        environmentCollisionGroup = 
            world->getConstraintSolver()->getCollisionDetector()->createCollisionGroup(env.get());

    } else {
        environmentCollisionGroup->addShapeFramesOf(env.get());
    }
}

void DartWorld::addObject(const SkeletonPtr& the_object){
    world->addSkeleton(the_object);
    object = the_object;
    objectCollisionGroup = world->getConstraintSolver()->getCollisionDetector()->
                createCollisionGroup(object.get());

}

void DartWorld::addRobot(DartManipulatorTemplate* the_robot){
    robot_ptr = the_robot;
    for (int i = 0; i < robot_ptr->bodies.size(); i++){
        world->addSkeleton(robot_ptr->bodies[i]);
    }
    
    robot_ptr->setupCollisionGroup(this->world);

    this->robotCollisionGroup = robot_ptr->mCollisionGroup;
}

void DartWorld::updateObjectPose(const Vector7d & object_pose){
    object->setPositions(pose7d_to_pose6d(object_pose));
}

void DartWorld::updateRobotConfig(const Eigen::VectorXd & robot_config){
    Vector7d object_pose = pose6d_to_pose7d(object->getPositions());
    robot_ptr->setConfig(robot_config, object_pose);
}

bool DartWorld::isRobotCollide(const Eigen::VectorXd & robot_config){
    updateRobotConfig(robot_config);
    dart::collision::CollisionOption option;
    dart::collision::CollisionResult result;
    bool collision = robotCollisionGroup->collide(environmentCollisionGroup.get(), option, &result);
    
    if (robot_ptr->ifCheckObjectCollision){
        bool object_collision = robotCollisionGroup->collide(objectCollisionGroup.get(), option, &result);
        return (object_collision | collision);
    }
    
    return collision;
}

void DartWorld::getObjectContacts(std::vector<ContactPoint>* pts){

    dart::collision::CollisionOption option;
    dart::collision::CollisionResult result;
    bool collision = objectCollisionGroup->collide(environmentCollisionGroup.get(), option, &result);


    Eigen::Vector6d object_pose = object->getPositions();
    Eigen::Matrix4d T;
    T = SE3Inv(pose2SE3(pose6d_to_pose7d(object_pose)));
    Eigen::Matrix3d R;
    R = T.block(0,0,3,3);
    Eigen::Vector3d p;
    p = T.block(0,3,3,1);

    // print out contact points
    // std::cout << "If in collision: " << result.isCollision() << std::endl;
    // int min_idx = 0;
    // double min_v = 10000;
    std::vector<ContactPoint> pts_;
    for (size_t i = 0; i < result.getNumContacts(); i++){
        Contact cp = result.getContact(i);
        ContactPoint pt;
        pt.p = R*cp.point + p;
        pt.n = R*cp.normal;
        pt.d = cp.penetrationDepth;
        pts->push_back(pt);

        // if (pt.p.sum() < min_v){
        //     min_v = pt.p.sum();
        //     min_idx = i;
        // }
    }

    // // TODO: delete this when ss_mode_enumeration fixed
    // for(int i = min_idx; i < pts_.size(); i++){
    //     pts->push_back(pts_[i]);
    // }
    // for (int i = 0; i < min_idx; i++){
    //     pts->push_back(pts_[i]);
    // }

}

void DartWorld::getObjectContacts(std::vector<ContactPoint>* pts, const Vector7d & object_pose){
    updateObjectPose(object_pose);
    getObjectContacts(pts);
}

void DartWorld::startWindow(int* pargc, char** argv){
    window->objectCollisionGroup = objectCollisionGroup;
    window->robotCollisionGroup = robotCollisionGroup;
    window->environmentCollisionGroup = environmentCollisionGroup;
    window->object = this->object;
    window->robot = this->robot_ptr;
    glutInit(pargc, argv);
    window->initWindow(800, 600, "Collisions");
    glutMainLoop();
}

void DartWorld::setPlaybackTrajectory(const std::vector<Vector7d>& object_traj, const std::vector<VectorXd>& robot_traj){

    if(this->window->object_positions.size() > 0){
        this->window->object_positions.clear();
    }

    if(this->window->robot_configs.size() > 0){
        this->window->robot_configs.clear();
    }

    for (auto x: object_traj){
        this->window->object_positions.push_back(pose7d_to_pose6d(x));
    }

    for(auto q: robot_traj){
        this->window->robot_configs.push_back(q);
    }
    this->window->play_mode = PLAY_BACK;

}

void DartWorld::setObjectTrajectory(const std::vector<Vector7d>& object_traj){

    if(this->window->object_positions.size() > 0){
        this->window->object_positions.clear();
    }

    if(this->window->robot_configs.size() > 0){
        this->window->robot_configs.clear();
    }

    for (auto x: object_traj){
        this->window->object_positions.push_back(pose7d_to_pose6d(x));
    }

    this->window->play_mode = SHOW_OBJECT_ALL;

}

void DartWorld::setSurfacePoints(const std::vector<ContactPoint>& pts){
    this->window->pts = pts;
    this->window->mDrawPoints = true;
}
