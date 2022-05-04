#ifndef MANIPULATORS_DARTURDFMANIPULATORTEMP
#define MANIPULATORS_DARTURDFMANIPULATORTEMP
    #include "DartURDFManipulatorTemplate.h"
#endif

class DartKR5: public virtual DartURDFManipulatorTemplate {

public:
    // config: [p1, n1] in the object frame
    DartKR5();
    DartKR5(const DartKR5& kr5): DartURDFManipulatorTemplate(kr5){}

    void setConfig(const VectorXd& config, const Vector7d& object_pose);

    void setupCollisionGroup(WorldPtr world);

    bool inverseKinematics(const VectorXd& mnp_config, const Vector7d& object_pose, VectorXd& result_joint_config);

    bool ifIKsolution(const VectorXd& mnp_config, const Vector7d& object_pose);

    void computeWorkspace(std::string save_to);

    void getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips);

    bool resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips);

    void Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, 
                                    std::vector<ContactPoint>* point_contacts);

    void applyJointPositions(const VectorXd& joint_config);

private:

    Vector6d reset_joint_config;


};