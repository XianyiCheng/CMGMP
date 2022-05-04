#include "DartPointManipulator.h"

#include "../utilities/sample.h"

#include "../contacts/contact_kinematics.h"



DartPointManipulator::DartPointManipulator(int n, double radius){
    
    this->n_pts = n;
    this->fingertip_radius = radius;

    for (int i=0; i< n; i++){
        SkeletonPtr ball = createFreeBall("ball_"+std::to_string(i), radius, Vector3d(0.3,0.3,0.8));
        this->addBody(ball);
    }
}

void DartPointManipulator::getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips){
    
    int n = int(config.size()/6);

    if (fingertips->size()!=0){
        fingertips->clear();
    }

    for (int i = 0; i < n; i++){
        fingertips->push_back(ContactPoint(config.segment<3>(6*i), config.segment<3>(6*i+3)));
    }

}

void DartPointManipulator::setConfig(const VectorXd& config, const Vector7d& object_pose){

    Eigen::Matrix4d T;
    T = pose2SE3(object_pose);
    Eigen::Matrix3d R;
    R = T.block(0,0,3,3);
    Eigen::Vector3d p;
    p = T.block(0,3,3,1);

    int n = int(config.size()/6);
    for (int i=0; i< this->n_pts; i++){
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        if (i < n){
            pos.tail(3) = R*Vector3d(config[6*i], config[6*i+1], config[6*i+2]) + p;
        }
        this->bodies[i]->setPositions(pos);
    }
}


bool DartPointManipulator::resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips){

    // n_on: number of relocating fingers 
    if (n_on == 0){
        return false;
    }
    
    std::vector<ContactPoint> fingertips;
    Vector7d x_o;
    getFingertipsOnObject(config, x_o, &fingertips);

    std::vector<ContactPoint> new_samples;
    // TODO: this can be again formulated as a search problem
    int n_pts = this->n_pts;
    int n_located = fingertips.size();
    int n_unlocated = n_pts - n_located;
    int N = object_surface.size();

    std::vector<ContactPoint> remain_mnps;
    std::vector<int> relocate_idxes;

    for(int i = 0; i < std::min(n_on, n_unlocated); i++){
        if (randd() > 0.5){
            // randomly choose to not locate this unlocated finger
            continue;
        }
        int idx = randi(N);
        new_samples.push_back(object_surface[idx]);
        remain_mnps.push_back(object_surface[idx]);
    }

    for(int i = 0; i < n_on - n_unlocated; i++){
        // randomly choose fingers to relocate

        int idx = randi(n_located);

        while(std::find(relocate_idxes.begin(), relocate_idxes.end(), idx) != relocate_idxes.end()){
            idx = randi(n_located);
        }

        // randomly choose to release the finger 
        if(randd() > 0.5){
            relocate_idxes.push_back(idx);
        } else {
            new_samples.push_back(object_surface[randi(N)]);
            relocate_idxes.push_back(idx);
        }


    }
    for (int k = 0; k < n_located; k++){
        // find remaining fingers
        if(std::find(relocate_idxes.begin(), relocate_idxes.end(), k) == relocate_idxes.end()){
            remain_mnps.push_back(fingertips[k]);
            new_samples.push_back(fingertips[k]);
        }
    }
    
    // copy_points(new_samples, new_mnps);
    copy_points(remain_mnps, remain_fingertips);

    // new configuration
    new_config.resize(new_samples.size()*6);
    for (int i = 0; i < new_samples.size(); i++){
        new_config.segment(6*i,3) = new_samples[i].p;
        new_config.segment(6*i+3,3) = new_samples[i].n;
    }

    new_samples.clear();
    return true;        

}

void DartPointManipulator::Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips,
 std::vector<ContactPoint>* point_contacts){
  
     copy_points(fingertips, point_contacts);
    
 }

void DartPointManipulator::setupCollisionGroup(WorldPtr world){
    auto collisionEngine = world->getConstraintSolver()->getCollisionDetector();
    this->mCollisionGroup = collisionEngine->createCollisionGroup();
    for (int i = 0; i < this->bodies.size(); i++){
        this->mCollisionGroup->addShapeFramesOf(this->bodies[i].get());
    }
}