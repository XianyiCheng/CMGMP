
#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "../contacts/contact_constraints.h"
#endif

#ifndef MANIPULATORS_MANIPULATORTEMPLATE   
#define MANIPULATORS_MANIPULATORTEMPLATE   
    #include "../manipulators/ManipulatorTemplate.h"
#endif

// wraps the simulation interface

// the world consists of: the environment one movable object, one manipulator

class WorldTemplate
{
public:
    
    WorldTemplate(){};
    WorldTemplate(const WorldTemplate& wt){};

    virtual void updateObjectPose(const Vector7d & object_pose) = 0;

    virtual void updateRobotConfig(const Eigen::VectorXd & robot_config) = 0;

    virtual bool isRobotCollide(const Eigen::VectorXd & robot_config) = 0;

    virtual void getObjectContacts(std::vector<ContactPoint>* pts) = 0;

    virtual void getObjectContacts(std::vector<ContactPoint>* pts, const Vector7d & object_pose) = 0;

    virtual void startWindow(int* pargc, char** argv) = 0;

    virtual void setPlaybackTrajectory(const std::vector<Vector7d>& object_traj, const std::vector<VectorXd>& robot_traj) = 0;

    virtual void setObjectTrajectory(const std::vector<Vector7d>& object_traj) = 0;

    virtual void setSurfacePoints(const std::vector<ContactPoint>& pts) = 0;

    virtual ManipulatorTemplate* getRobot() = 0;

};