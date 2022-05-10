#include "DartKR5.h"
#include "../utilities/sample_grasp.h"
#include "../contacts/contact_kinematics.h"

Matrix3d znormal2rotmat(Vector3d normal){
    Vector3d z(0,0,1); 
    Vector3d rotvec = z.cross(normal);
    Matrix3d R;

    if ((rotvec.norm() < 1e-6) && (signbit(z[2]) == signbit(normal[2]))){    
        R.setIdentity();

    } else if ((rotvec.norm() < 1e-6) && (signbit(z[2]) != signbit(normal[2]))){
        R << 1,0,0,0,-1,0,0,0,1;

    }
    else{
        rotvec = rotvec*(1/rotvec.norm());
        double rotangle = acos(z.dot(normal)/normal.norm());
        AngleAxisd aaxis(rotangle, rotvec);
        R = aaxis.toRotationMatrix();
    }    
    return R;
}

DartKR5::DartKR5(){
    this->loadManipulator("KR5", "/home/xianyi/Research/CMGMP/data/urdf/KR5/KR5 sixx R650.urdf");
    this->n_pts = 1;
    this->ifCheckObjectCollision = true;
    this->NumDofs = this->bodies[0]->getNumJoints();

    this->reset_joint_config << 0, -0.01, 0.24, 0, 1.2 , 0;
    std::cout << "Robot body nodes: " << std::endl;
    for (int i = 0; i < this->bodies[0]->getNumBodyNodes(); i++){
        BodyNode* bn = this->bodies[0]->getBodyNode(i);
        std::cout << bn->getName() << ", ";
    }
    std::cout << std::endl;

    // add a point manipulator

    double e_r = 0.005;
    double e_h = 0.07;

    {

        WeldJoint::Properties properties;
        properties.mName = "ee_joint";
        properties.mT_ParentBodyToJoint = Eigen::Isometry3d::Identity();
        properties.mT_ChildBodyToJoint = Eigen::Isometry3d::Identity();
        properties.mT_ChildBodyToJoint.linear() = (AngleAxisd(0.5*M_PI, Vector3d::UnitY())).matrix();
        properties.mT_ChildBodyToJoint.translation() = Vector3d(0,0,0.5*e_h);


        BodyNode* parent = this->bodies[0]->getBodyNode(this->bodies[0]->getNumBodyNodes()-1);
        BodyNode* bn = this->bodies[0]->createJointAndBodyNodePair<WeldJoint>(
            parent, properties, BodyNode::AspectProperties("end_effector")).second;

        std::shared_ptr<CylinderShape> shape =std::make_shared<CylinderShape>(e_r, e_h);
        auto shapeNode
            = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
                shape);
        shapeNode->getVisualAspect()->setColor(Eigen::Vector3d(0.3,0.5,0.3));
        shapeNode->getVisualAspect()->setAlpha(0.85);
    }
    

    {
        WeldJoint::Properties properties;
        properties.mName = "ee_balljoint";
        properties.mT_ParentBodyToJoint = Eigen::Isometry3d::Identity();
        properties.mT_ChildBodyToJoint = Eigen::Isometry3d::Identity();
        properties.mT_ChildBodyToJoint.translation() = Vector3d(0,0,0.5*e_h + e_r);


        BodyNode* parent = this->bodies[0]->getBodyNode(this->bodies[0]->getNumBodyNodes()-1);
        BodyNode* bn = this->bodies[0]->createJointAndBodyNodePair<WeldJoint>(
            parent, properties, BodyNode::AspectProperties("end_effector_ball")).second;

        std::shared_ptr<SphereShape> shape =std::make_shared<SphereShape>(e_r);
        auto shapeNode
            = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
                shape);
        shapeNode->getVisualAspect()->setColor(Eigen::Vector3d(0.3,0.5,0.3));
        shapeNode->getVisualAspect()->setAlpha(0.85);

    }

    // set IK
    this->ik = this->bodies[0]->getBodyNode("end_effector_ball")->getIK(true);

    // this->computeWorkspace("./KR5_workspace.csv");

    Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
    tf.translation() = Eigen::Vector3d(0, 0.0, 0.0);
    this->bodies[0]->getJoint(0)->setTransformFromParentBodyNode(tf);

}

void DartKR5::applyJointPositions(const VectorXd& joint_config){
    this->bodies[0]->setPositions(joint_config);
}

void DartKR5::setConfig(const VectorXd& config, const Vector7d& object_pose){
    if(config.size() < 6){
        this->bodies[0]->setPositions(this->reset_joint_config);
        return;
    }
    VectorXd joint_config;
    bool ifIK = inverseKinematics(config, object_pose, joint_config);
    if (ifIK){
        this->bodies[0]->setPositions(joint_config);
    } else {
        printf("Error! Cannot find any IK solution.\n");
    }
}

void DartKR5::setupCollisionGroup(WorldPtr world){
    auto collisionEngine = world->getConstraintSolver()->getCollisionDetector();
    this->mCollisionGroup = collisionEngine->createCollisionGroup();

    // not considering the collision of "base_link" and "end_effector_ball"
    
    for (int i = 1; i < this->bodies[0]->getNumBodyNodes() - 1; i++){
        this->mCollisionGroup->addShapeFramesOf(
            this->bodies[0]->getBodyNode(i));
    }

}

bool DartKR5::inverseKinematics(const VectorXd& config, const Vector7d& object_pose, VectorXd& result_joint_config){

    if(config.size() < 6){
        return true;
    }


    Eigen::Matrix4d T;
    T = pose2SE3(object_pose);
    Eigen::Matrix3d R;
    R = T.block(0,0,3,3);
    Eigen::Vector3d p;
    p = T.block(0,3,3,1);

    Vector6d config_world;
    config_world.head(3) = R*config.head(3) + p;
    config_world.tail(3) = R*config.tail(3);

    auto targetFrame
      = SimpleFrame::createShared(Frame::World());
    targetFrame->setRotation(znormal2rotmat(-config_world.tail(3)));
    // targetFrame->setRotation(znormal2rotmat(Vector3d(-0.5,0, 0.866025)));
    targetFrame->setRotation(znormal2rotmat(Vector3d(-0.17364817766693033,0, 0.984807753012208)));
    // targetFrame->setRotation(znormal2rotmat(Vector3d(0,0,1)));


    targetFrame->setTranslation(config_world.head(3));
    

    this->ik->setTarget(targetFrame);
    bool ifIKfound = this->ik->findSolution(result_joint_config);

    // Eigen::Matrix3d COM = this->bodies[0]->getBodyNode("end_effector_ball")->getTransform().linear();// * Vector3d(0,0,0);
    // std::cout << COM << std::endl;
    return ifIKfound;
}

bool DartKR5::ifIKsolution(const VectorXd& mnp_config, const Vector7d& object_pose){
    VectorXd jc;
    return this->inverseKinematics(mnp_config, object_pose, jc);
}

void DartKR5::computeWorkspace(std::string save_to){
    
    // draw points

    MatrixXd data(5002, 3);

    VectorXd lb = this->bodies[0]->getPositionLowerLimits();
    VectorXd ub = this->bodies[0]->getPositionUpperLimits();

    this->bodies[0]->setPositions(lb);
    {
        Vector3d COM = this->bodies[0]->getBodyNode("end_effector_ball")->getTransform().translation();// * Vector3d(0,0,0);
        data.row(0)= COM.transpose();
        // show pt
        SkeletonPtr ball = createFreeBall("target" + std::to_string(0), 0.005, Vector3d(0.3,0.3,0.8));
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        pos.tail(3) = COM;
        ball->setPositions(pos);
        this->addBody(ball);
    }

    this->bodies[0]->setPositions(ub);
    {
        Vector3d COM = this->bodies[0]->getBodyNode("end_effector_ball")->getTransform().translation();// * Vector3d(0,0,0);
        data.row(1) = COM.transpose();

        // show pt
        SkeletonPtr ball = createFreeBall("target" + std::to_string(1), 0.005, Vector3d(0.3,0.3,0.8));
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        pos.tail(3) = COM;
        ball->setPositions(pos);
        this->addBody(ball);
    }

    for (int i = 2; i < 5000; i++){
        for (int k = 0; k < 6; k++){
            this->bodies[0]->setPosition(k, lb[k] + randd()*(ub[k] - lb[k]));
        }
        Vector3d COM = this->bodies[0]->getBodyNode("end_effector_ball")->getTransform().translation();// * Vector3d(0,0,0);
        data.row(i) = COM.transpose();

        // show pt
        SkeletonPtr ball = createFreeBall("target" + std::to_string(i), 0.005, Vector3d(0.3,0.3,0.8));
        Eigen::Vector6d pos(Eigen::Vector6d::Zero());
        pos.tail(3) = COM;
        ball->setPositions(pos);
        this->addBody(ball);
    }
    saveData(save_to, data);
}

void DartKR5::getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips)
{
    if (config.size() < 6){
        if (fingertips->size() > 0){
            fingertips->clear();
        }
        return;
    }

    ContactPoint pt;
    pt.p = config.head(3);
    pt.n = config.tail(3);
    pt.d = 0;
    fingertips->push_back(pt);

}

bool DartKR5::resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, 
                                const std::vector<ContactPoint>& object_surface,
                                VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips)
{
    // n_on: number of relocating fingers 
    if (n_on == 0 || config.size()!=0){
        // if (config.size() < 6){
        //     return false;
        // }
        // Matrix3d R;
        // double angle;
        // Vector3d axis;
        
        // for (int k = 0; k < 4; k++)
        // {
        //     if (k == 1 || k == 3){
        //         angle = -M_PI*15/180;
        //     } else {
        //         angle = M_PI*15/180;
        //     }
        //     if (k == 2 || k == 3){
        //         axis = Eigen::Vector3d::UnitY();
        //     } else {
        //         axis = Eigen::Vector3d::UnitX();
        //     }

        //     AngleAxisd aaxis(angle, axis);
        //     R = aaxis.toRotationMatrix();

        //     new_config.head(3) = config.head(3);
        //     new_config.tail(3) = R*config.tail(3);

        //     bool isIK = this->ifIKsolution(new_config, object_pose);

        //     if (isIK){
        //         return true;
        //     } 
        // }
        new_config.resize(0);
        return true;
    }
    
    if (remain_fingertips->size() > 0){
        remain_fingertips->clear();
    }

    int n_sample = 30;
    int i_sample = 0;
    int N_sample = object_surface.size();

    while (i_sample < n_sample){
        i_sample++;
        int k = randi(N_sample);
        VectorXd cf(6);

        // AngleAxisd aaxis((randd()-0.5)*0.5*M_PI, 
        //                 AngleAxisd(randd()*2*M_PI, Eigen::Vector3d::UnitZ()).matrix()*Eigen::Vector3d::UnitX());
        // Matrix3d R = aaxis.toRotationMatrix();

        cf.head(3) = object_surface[k].p;
        cf.tail(3) = object_surface[k].n;
        
        bool isIK = this->ifIKsolution(cf, object_pose);

        if (isIK){
            new_config = cf;
            return true;
        }
    }
    printf("Failed to sample a manipulator contact! Please check your IK. \n");

    return false;
 
}


void DartKR5::Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, 
                                std::vector<ContactPoint>* point_contacts)
{
    double r = 0.005;

    Matrix3d xs;
    xs << 1, 0, 0,
     -0.5, 0.865, 0,
     -0.5, -0.865, 0;

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
    return;
}



