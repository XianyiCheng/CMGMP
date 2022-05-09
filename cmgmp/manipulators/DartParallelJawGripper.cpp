#include "DartParallelJawGripper.h"

#include "../utilities/sample.h"

#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "../contacts/contact_kinematics.h"
#endif

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

Quaterniond xnormal2quat(Vector3d normal){
    Vector3d z(1,0,0); 
    Vector3d rotvec = z.cross(normal);
    Matrix3d R;

    Quaterniond q;

    if ((rotvec.norm() < 1e-6) && (signbit(z[0]) == signbit(normal[0]))){    
        R.setIdentity();
        q = Quaterniond(R);
    } else if ((rotvec.norm() < 1e-6) && (signbit(z[0]) != signbit(normal[0]))){
        R << -1,0,0,0,1,0,0,0,-1;
        q = Quaterniond(R);
    }
    else{
        rotvec = rotvec*(1/rotvec.norm());
        double rotangle = acos(z.dot(normal)/normal.norm());
        AngleAxisd aaxis(rotangle, rotvec);
        // R = aaxis.toRotationMatrix();
        q = Quaterniond(aaxis);
    }    
    return q;
}

Quaterniond ynormal2quat(Vector3d normal){
    Vector3d z(0,1,0); 
    Vector3d rotvec = z.cross(normal);
    Matrix3d R;

    Quaterniond q;

    if ((rotvec.norm() < 1e-6) && (signbit(z[1]) == signbit(normal[1]))){    
        R.setIdentity();
        q = Quaterniond(R);
    } else if ((rotvec.norm() < 1e-6) && (signbit(z[1]) != signbit(normal[1]))){
        R << 1,0,0,0,-1,0,0,0,-1;
        q = Quaterniond(R);
    }
    else{
        rotvec = rotvec*(1/rotvec.norm());
        double rotangle = acos(z.dot(normal)/normal.norm());
        AngleAxisd aaxis(rotangle, rotvec);
        // R = aaxis.toRotationMatrix();
        q = Quaterniond(aaxis);
    }    
    return q;
}

Quaterniond znormal2quat(Vector3d normal){
    Vector3d z(0,0,1); 
    Vector3d rotvec = z.cross(normal);
    Matrix3d R;

    Quaterniond q;

    if ((rotvec.norm() < 1e-6) && (signbit(z[2]) == signbit(normal[2]))){    
        R.setIdentity();
        q = Quaterniond(R);
    } else if ((rotvec.norm() < 1e-6) && (signbit(z[2]) != signbit(normal[2]))){
        R << 1,0,0,0,-1,0,0,0,1;
        q = Quaterniond(R);
    }
    else{
        rotvec = rotvec*(1/rotvec.norm());
        double rotangle = acos(z.dot(normal)/normal.norm());
        AngleAxisd aaxis(rotangle, rotvec);
        // R = aaxis.toRotationMatrix();
        q = Quaterniond(aaxis);
    }    
    return q;
}

Vector3d sample_unit_vector(){
    Vector3d ub(1,1,1);
    Vector3d lb(-1,-1,-1);
    Vector3d p;
    p = sample_position(ub, lb);
    p = p/p.norm();
    return p;
}

Vector3d sample_orthogonal_unit_vector(const Vector3d& axis){
    Vector3d ub(1,1,1);
    Vector3d lb(-1,-1,-1);
    Vector3d p;
    p = axis.cross(sample_position(ub, lb));
    p = p/p.norm();
    return p;
}

DartParallelJawGripper::DartParallelJawGripper(double max_d, double radius, double fl){
    
    this->ifCheckObjectCollision = true;
    this->max_open_dist = max_d;
    this->finger_length = fl;
    this->fingertip_radius = radius;

    if (this->finger_length == 0){
        this->finger_length = 0.6*this->max_open_dist;
    }


    this->n_pts = 2;

    // add fingertips
    for (int i=0; i< this->n_pts; i++){
        SkeletonPtr ball = createFreeBall("fingertip_"+std::to_string(i), radius, Vector3d(0.3,0.3,0.8));
        this->addBody(ball);
    }

    // add fingers
    for (int i=0; i< this->n_pts; i++){
        SkeletonPtr cylinder = createFreeCylinder("finger_"+std::to_string(i), radius*0.7, this->finger_length - 2.2*radius);
        this->addBody(cylinder);
    }
    // add jaw
    SkeletonPtr cylinder = createFreeCylinder("jaw", radius*0.8, this->max_open_dist);
    this->addBody(cylinder);
}

bool DartParallelJawGripper::setupCollisionGroup(WorldPtr world){
    
    // don't collide the fingertips.
    auto collisionEngine = world->getConstraintSolver()->getCollisionDetector();
    this->mCollisionGroup = collisionEngine->createCollisionGroup(
        this->bodies[2].get(), this->bodies[3].get(), this->bodies[4].get());
}

bool DartParallelJawGripper::ifIK(const VectorXd& config){
    if (config.size() < 15){
        return true;
    } else {

        Vector3d p1(config[0], config[1], config[2]);
        Vector3d p2(config[6], config[7], config[8]);

        double d = (p1 - p2).norm();

        if (d > this->max_open_dist || d < 1e-4){
            return false;
        } else {
            return true;
        }
    }
}

bool DartParallelJawGripper::ifIK(const std::vector<ContactPoint>& fingertips){
    if (fingertips.size()<2){
        return true;
    }

    double d = (fingertips[0].p - fingertips[1].p).norm();

    if (d > this->max_open_dist || d < 1e-4){
        return false;
    } else {
        return true;
    }
}

void DartParallelJawGripper::getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips){
    
    if (config.size() < 3){
        if (fingertips->size()!=0){
        fingertips->clear();
        }
        return;
    }

    int n = int((config.size()-3)/6);

    if (fingertips->size()!=0){
        fingertips->clear();
    }

    for (int i = 0; i < n; i++){
        fingertips->push_back(ContactPoint(config.segment<3>(6*i), config.segment<3>(6*i+3)));
    }

}

void DartParallelJawGripper::setConfig(const VectorXd& config, const Vector7d& object_pose){

    int n = 0;
    double r = 1.1*this->fingertip_radius;
    // if there is no config
    if (config.size() < 3){
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        pos.tail(3) = Vector3d(100,100,100);
        for(int i=0; i < this->bodies.size(); i++){
            this->bodies[i]->setPositions(pos);
        }
        return;
    } else {
        n = int((config.size()-3)/6);
    }


    Eigen::Matrix4d T;
    T = pose2SE3(object_pose);
    Eigen::Matrix3d R;
    R = T.block(0,0,3,3);
    Eigen::Vector3d p;
    p = T.block(0,3,3,1);

    // set fingertips locations
    for (int i=0; i< this->n_pts; i++){
        Vector3d pos;
        if (n==1){
            pos = Vector3d(config[0] - r*config[3], config[1] - r*config[4], config[2] - r*config[5]);
        } else {
            pos = Vector3d(config[6*i] - r*config[6*i+3], config[6*i+1]- r*config[6*i+4], config[6*i+2]- r*config[6*i+5]);
        }
        Eigen::Vector6d pos_w(Eigen::Vector6d::Zero());
        
        pos_w.tail(3) = R*pos + p;
        
        this->bodies[i]->setPositions(pos_w);
    }

    // set finger positions
    Vector3d gripper_pos;
    Vector3d finger_pos_1;
    Vector3d finger_pos_2;
    Vector3d nt; // gripper normal
    Vector3d nf; // finger normal
    
    if (n == 1){
        nt = Vector3d(config[6], config[7], config[8]);
        Vector3d nj(config[9], config[10], config[11]);
        gripper_pos = Vector3d(config[0] - r*config[3], config[1] - r*config[4], config[2] - r*config[5]); - this->finger_length*nj;
        finger_pos_1 = Vector3d(config[0] - r*config[3], config[1] - r*config[4], config[2] - r*config[5]); - 0.5*this->finger_length*nj;
        finger_pos_2 = finger_pos_1;
        nf << nj[0], nj[1], nj[2];
    } else if (n == 2){
        nt = Vector3d(config[0]-config[6], config[1]-config[7], config[2]-config[8]);
        Vector3d nj(config[12], config[13], config[14]);
        gripper_pos = Vector3d((config[0]+config[6])/2, 
            (config[1]+config[7])/2, (config[2]+config[8])/2) - this->finger_length*nj;
        finger_pos_1 = Vector3d(config[0] - r*config[3], config[1] - r*config[4], config[2] - r*config[5]) - 0.5*this->finger_length*nj;
        finger_pos_2 = Vector3d(config[6] - r*config[9], config[7] - r*config[10], config[8] - r*config[11]) - 0.5*this->finger_length*nj;
        nf << nj[0], nj[1], nj[2];
    }

    // jaw
    Quaterniond qt = znormal2quat(nt);
    {
        Vector7d gripper_pose;
        gripper_pose << gripper_pos[0], gripper_pos[1], gripper_pos[2], qt.x(), qt.y(), qt.z(), qt.w();
        Matrix4d gripper_T = T*pose2SE3(gripper_pose);
        this->bodies[4]->setPositions(pose7d_to_pose6d(SE32pose(gripper_T)));
    }

    // fingers
    Quaterniond qf = znormal2quat(nf);
    {
        Vector7d finger_pose;
        finger_pose << finger_pos_1[0], finger_pos_1[1], finger_pos_1[2], qf.x(), qf.y(), qf.z(), qf.w();
        Matrix4d finger_T = T*pose2SE3(finger_pose);
        this->bodies[2]->setPositions(pose7d_to_pose6d(SE32pose(finger_T)));
    }
    {
        Vector7d finger_pose;
        finger_pose << finger_pos_2[0], finger_pos_2[1], finger_pos_2[2], qf.x(), qf.y(), qf.z(), qf.w();
        Matrix4d finger_T = T*pose2SE3(finger_pose);
        this->bodies[3]->setPositions(pose7d_to_pose6d(SE32pose(finger_T)));
    }


}

bool DartParallelJawGripper::resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips)

{

    // n_on: number of relocating fingers 
    
    
   std::vector<ContactPoint> new_samples;
   std::vector<ContactPoint> fingertips;

   Vector7d xo;

   this->getFingertipsOnObject(config, xo, &fingertips);

    // TODO: this can be again formulated as a search problem
    int n_pts = this->n_pts;
    int n_located = fingertips.size();
    int n_unlocated = n_pts - n_located;
    int N = object_surface.size();

    if (n_on == 0 && n_located!=0){
        if (n_located == 1){
            new_config.resize(12);
            new_config.segment(0,6) = config.segment(0,6);
            // sample jaw direction
            Vector3d n_t = sample_unit_vector();
            // n_t = -new_samples[0].n;
            // sample fingertip direction
            Vector3d n_f = sample_orthogonal_unit_vector(n_t);
            new_config.segment(6,3) = n_t;
            new_config.segment(9,3) = n_f;
            copy_points(fingertips, remain_fingertips);
        } else { // n_located == 2
            new_config.resize(15);
            Vector3d n_t = fingertips[0].p - fingertips[1].p;
            // sample fingertip direction
            Vector3d n_f = sample_orthogonal_unit_vector(n_t);
            new_config.segment(0,12) = config.segment(0,12);
            new_config.segment(12,3) = n_f;
            copy_points(fingertips, remain_fingertips);
        }
        return true;
    }

    if (n_on == 0 && n_located==0){
        n_on = randi(n_pts)+1;
    }

    // n_on = 2;

    int samples = 10;
    
    for(int k_sample = 0; k_sample < samples; k_sample++){
        
        std::vector<ContactPoint> remain_mnps;
        std::vector<int> relocate_idxes;

        for(int i = 0; i < std::min(n_on, n_unlocated); i++){
            int idx = randi(N);
            new_samples.push_back(object_surface[idx]);
            remain_mnps.push_back(object_surface[idx]);
        }

        // check IK during a special case of relocation
        if (n_on==2 && n_unlocated==1){
            std::vector<ContactPoint> temp_mnps;
            temp_mnps.push_back(remain_mnps[0]);
            temp_mnps.push_back(fingertips[0]);
            if (!this->ifIK(temp_mnps)){
                new_samples.clear();
                remain_mnps.clear();
                continue;
            }
        }

        for(int i = 0; i < n_on - n_unlocated; i++){
            // randomly choose fingers to relocate
            int idx = randi(n_located);
            while(std::find(relocate_idxes.begin(), relocate_idxes.end(), idx) != relocate_idxes.end()){
                idx = randi(n_located);
            }
            new_samples.push_back(object_surface[randi(N)]);
            relocate_idxes.push_back(idx);
        }
        for (int k = 0; k < n_located; k++){
            if(std::find(relocate_idxes.begin(), relocate_idxes.end(), k) == relocate_idxes.end()){
                remain_mnps.push_back(fingertips[k]);
                new_samples.push_back(fingertips[k]);
            }
        }

        
        bool isIK = this->ifIK(new_samples);

        if (isIK){
            Vector3d n_t;
            Vector3d n_f;
            if (new_samples.size()==1){
                new_config.resize(12);
                new_config.segment(0,3) = new_samples[0].p;
                new_config.segment(3,3) = new_samples[0].n;
                // sample jaw direction
                n_t = sample_unit_vector();
                // n_t = -new_samples[0].n;
                // sample fingertip direction
                n_f = sample_orthogonal_unit_vector(n_t);
                new_config.segment(6,3) = n_t;
                new_config.segment(9,3) = n_f;
            } else if(new_samples.size()==2){
                new_config.resize(15);
                n_t = new_samples[0].p - new_samples[1].p;
                // sample fingertip direction
                n_f = sample_orthogonal_unit_vector(n_t);
                new_config.segment(0,3) = new_samples[0].p;
                new_config.segment(3,3) = new_samples[0].n;
                new_config.segment(6,3) = new_samples[1].p;
                new_config.segment(9,3) = new_samples[1].n;
                new_config.segment(12,3) = n_f;
            }

            copy_points(remain_mnps, remain_fingertips);
            new_samples.clear();
            remain_mnps.clear();
            return true;
        } 
        new_samples.clear();
        remain_mnps.clear();
    }
    
    return false; 

} 

void DartParallelJawGripper::Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, 
                                std::vector<ContactPoint>* point_contacts)
{
    Matrix3d xs;
    xs << 1, 0, 0,
     -0.5, 0.865, 0,
     -0.5, -0.865, 0;

    double r = this->fingertip_radius;
    for(auto& pt: fingertips){
        Matrix6d Adgco = contact_jacobian(pt.p, pt.n);
        Matrix3d Roc = Adgco.topLeftCorner(3, 3).transpose();
        for (int i = 0; i < 3; i++){
            ContactPoint gpt;
            gpt.p = pt.p + r*Roc*(xs.row(i)).transpose();
            gpt.n = pt.n;
            gpt.d = 0;
            point_contacts->push_back(gpt);
        }
    }
}