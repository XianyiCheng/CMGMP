#include "DartDDHandFreeY.h"

#include "../utilities/sample.h"

#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "../contacts/contact_kinematics.h"
#endif

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

void DartDDHandFreeY::setConfig(const VectorXd& config, const Vector7d& object_pose){


    double y = 0.0;

    // first set all fingers to its reset position
    Eigen::Vector6d pos(Eigen::Vector6d::Zero());
    pos.tail(3) = hand_pose.head(3);
    for (int i = 0; i < this->bodies.size(); i++){
        this->bodies[i]->setPositions(pos);
    }

    this->obj_pose = object_pose;

    Matrix3d R_WH = quat2rotm(hand_pose[3], hand_pose[4], hand_pose[5], hand_pose[6]);

    Matrix4d T = pose2SE3(object_pose);
    Matrix3d R_WO = T.block(0,0,3,3);
    Vector3d p_WO = T.block(0,3,3,1);

    Vector3d finger_axis;

    if (fingertype == I_FINGER){
        finger_axis = R_WH*Vector3d(0,0,1);
    } else if (fingertype == L_FINGER) {
        finger_axis = R_WH*Vector3d(0,1,0);
    }

    Vector3d orn = znormal2aaxis(finger_axis);

    {
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        pos.head(3) = orn;

        std::vector<int> idx;
        getActiveFingerIdx(config, &idx);

        for (auto i:idx) {
            if (fingertype == I_FINGER){

                Vector3d pp = R_WO*Vector3d(config[6*i+0], config[6*i+1], config[6*i+2]) + p_WO;

                y = pp[1];

                pos.tail(3) = pp + radius*1.1*finger_axis;
                this->bodies[i]->setPositions(pos);
                pos.tail(3) = pp + 0.5*(finger_length + 1.2*radius)*finger_axis;
                this->bodies[i+2]->setPositions(pos);

            } else if (fingertype == L_FINGER) {
                Vector3d pp = R_WO*(Vector3d(config[6*i+0], config[6*i+1], config[6*i+2]) 
                    - radius*Vector3d(config[6*i+3], config[6*i+4], config[6*i+5])) + p_WO;
                    
                pos.tail(3) = pp;
                this->bodies[i]->setPositions(pos);
            }
        }
    }

    this->hand_pose[1] = y;

    this->setHandFrameTransform(this->hand_pose);
}

VectorXd DartDDHandFreeY::getHandFrameConfig(const VectorXd& config, const Vector7d& object_pose){
    
    std::vector<int> idx;
    getActiveFingerIdx(config, &idx);

    std::vector<Vector2d> pts;
    

    pts.push_back(Vector2d(-0.04, -0.05));
    pts.push_back(Vector2d(0.04, -0.05));

    double y = 0.0;

    if (idx.size() > 0) {

        // compute p1 p2 in the hand frame
        Matrix4d T_WH = pose2SE3(hand_pose);
        Matrix4d T_WO = pose2SE3(object_pose);
        Matrix4d T_HO = SE3Inv(T_WH)*T_WO;

        Matrix3d R_HO = T_HO.block(0,0,3,3);
        Vector3d p_HO = T_HO.block(0,3,3,1);

        Matrix3d R_WO = T_WO.block(0,0,3,3);
        Vector3d p_WO = T_WO.block(0,3,3,1);        

        for (auto i: idx){
            // Vector3d pp = R_HO*(config.segment<3>(6*i) + 0.001*config.segment<3>(6*i+3)) + p_HO;
            Vector3d pp = R_HO*config.segment<3>(6*i) + p_HO;
            Vector3d pp_w = R_WO*config.segment<3>(6*1) + p_WO;

            y = pp_w[1];
            
            if (i == 0){

                pp[0] = pp[0] - dx;
                pp[2] = pp[2] + dz;

                pts[i] = Vector2d(pp[0], pp[2]);

            } else {

                pp[0] = pp[0] + dx;
                pp[2] = pp[2] + dz;

                pts[i] = Vector2d(pp[0], pp[2]);

            }

        }
    }

    VectorXd rc(5);
    rc << pts[0][0], pts[0][1], pts[1][0], pts[1][1], y;

    return rc;
}

bool DartDDHandFreeY::resampleFingers(int n_on, const VectorXd& config, 
    const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips){

    // n_on: number of relocating fingers 
    if (n_on == 0){
        return false;
    }

    int N = object_surface.size();

    std::vector<int> located_idx;
    getActiveFingerIdx(config, &located_idx);

    int n_located = located_idx.size();

    for (int sample = 0; sample < 20; sample++){
        
        std::vector<ContactPoint> remain_mnps;

        VectorXd new_sample(12);
        if (config.size() < 12){
            new_sample << 0,0,0,0,0,0,0,0,0,0,0,0;
        } else {
            new_sample = config;
        }

        if (n_on == 1){
            // randomly choose a finger to relocate
            int finger_idx = randi(2);
            
            if (n_located = 0){
                int sample_idx = randi(N);
                new_sample.segment<3>(6*finger_idx) = object_surface[sample_idx].p;
                new_sample.segment<3>(6*finger_idx+3) = object_surface[sample_idx].n;
                remain_mnps.push_back(object_surface[sample_idx]);
            } 
            else if (n_located == 1){

                if (finger_idx == located_idx[0]){
                    int sample_idx = randi(N);
                    new_sample.segment<3>(6*finger_idx) = object_surface[sample_idx].p;
                    new_sample.segment<3>(6*finger_idx+3) = object_surface[sample_idx].n;
                } else {
                    int sample_1 = 0;
                    while (sample_1 < 20){
                        sample_1 ++;
                        int sample_idx = randi(N);
                        if (abs(config[6*located_idx[0]+1] - object_surface[sample_idx].p[1]) < max_dy){
                            new_sample.segment<3>(6*finger_idx) = object_surface[sample_idx].p;
                            new_sample.segment<3>(6*finger_idx+3) = object_surface[sample_idx].n;
                            new_sample[6*finger_idx+1] = config[6*located_idx[0]+1];
                            remain_mnps.push_back(ContactPoint(config.segment<3>(6*located_idx[0]), config.segment<3>(6*located_idx[0]+3)));
                            break;
                        }
                    }
                    if (sample_1 >= 20){
                        continue;
                    }
                }

            } 
            else if (n_located == 2){

                int left_idx = 1 - finger_idx;
                
                int sample_1 = 0;

                while (sample_1 < 20){

                    sample_1 ++;
                    int sample_idx = randi(N);

                    if (abs(config[6*left_idx+1] - object_surface[sample_idx].p[1]) < max_dy){

                        new_sample.segment<6>(6*left_idx) = config.segment<6>(6*left_idx);

                        new_sample.segment<3>(6*finger_idx) = object_surface[sample_idx].p;
                        new_sample.segment<3>(6*finger_idx+3) = object_surface[sample_idx].n;

                        new_sample[6*finger_idx+1] = config[6*left_idx+1];

                        remain_mnps.push_back(ContactPoint(config.segment<3>(6*left_idx), config.segment<3>(6*left_idx+3)));

                        break;
                    }
                }
                if (sample_1 >= 20){
                    continue;
                }

            }
        }

        if (n_on == 2){

            int sample_1 = 0;

            while (sample_1 < 20){
            
                int sample_idx = randi(N);
                int sample_idx1 = randi(N);

                while(sample_idx1 == sample_idx && N > 1){
                    sample_idx1 = randi(N);
                }

                if (abs(object_surface[sample_idx].p[1] - object_surface[sample_idx1].p[1]) < max_dy){

                    new_sample.segment<3>(0) = object_surface[sample_idx].p;
                    new_sample.segment<3>(3) = object_surface[sample_idx].n;

                    new_sample.segment<3>(6) = object_surface[sample_idx1].p;
                    new_sample.segment<3>(9) = object_surface[sample_idx1].n;

                    new_sample[1] = (new_sample[1] + new_sample[7])/2;
                    new_sample[7] = new_sample[1];

                    // if the finger is not located before
                
                    if (std::find(located_idx.begin(), located_idx.end(), 0) == located_idx.end()){
                        remain_mnps.push_back(object_surface[sample_idx]);
                    }

                    if (std::find(located_idx.begin(), located_idx.end(), 1) == located_idx.end()){
                        remain_mnps.push_back(object_surface[sample_idx1]);
                    }

                    break;

                }
            }
        
            if (sample_1 >= 20){
                continue;
            }
        }


        if(this->isIK(new_sample, object_pose)){
            new_config = new_sample;
            copy_points(remain_mnps, remain_fingertips);
            return true;
        }
    }

    printf("Failed to sample a manipulator contact! Please check your IK. \n");
    return false;    

}