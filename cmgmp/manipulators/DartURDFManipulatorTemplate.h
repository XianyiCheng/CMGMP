#ifndef MANIPULATORS_DARTMANIPULATORTEMPLATE
#define MANIPULATORS_DARTMANIPULATORTEMPLATE
    #include "DartManipulatorTemplate.h"
#endif

#ifndef DART_UTILS
#define DART_UTILS
    #include "../dart_utils/dart_utils.h"
#endif

// template for loading manipulators from urdf files

// set contact points
// inverse kinematics
// forward kinematics

//

class DartURDFManipulatorTemplate: public virtual DartManipulatorTemplate {

public:
    // config should be defined wrt to the object
    int NumDofs = 0;
    InverseKinematicsPtr ik = 0;
    DartURDFManipulatorTemplate(){}
    DartURDFManipulatorTemplate(const DartURDFManipulatorTemplate& dumt): DartManipulatorTemplate(dumt){}

    void loadManipulator(const std::string & name, const std::string & filePath){
        SkeletonPtr robot = createRobot(name, filePath);
        for(size_t i=0; i < robot->getNumJoints(); i++){
            robot->getJoint(i)->setPositionLimitEnforced(true);
        }
        robot->enableSelfCollisionCheck();
        this->addBody(robot);
    }
    
    virtual bool inverseKinematics(const VectorXd& mnp_config, const Vector7d& object_pose, VectorXd& result_joint_config) = 0;

    virtual void computeWorkspace(std::string save_to) = 0;

};