#ifndef MANIPULATORS_DARTMANIPULATORTEMPLATE
#define MANIPULATORS_DARTMANIPULATORTEMPLATE
    #include "DartManipulatorTemplate.h"
#endif

#ifndef DART_UTILS
#define DART_UTILS
    #include "../dart_utils/dart_utils.h"
#endif

#include "ddhand/handik.h"

#define L_FINGER 0
#define I_FINGER 1

#define FINGER_0 0
#define FINGER_1 1

Eigen::Vector3d quat2aaxis(double x, double y, double z, double w);

Matrix3d quat2rotm(double x, double y, double z, double w);

Vector3d znormal2aaxis(Vector3d normal);

void getActiveFingerIdx(const VectorXd& config, std::vector<int>* idx);

class DartDDHand: public virtual DartManipulatorTemplate {

public:

    int fingertype = -1;
    double radius = 0.005;
    double finger_length = 0.02;

    double dz = 28.63*1e-3;
    double dx = 15.64*1e-3;

    Vector7d hand_pose;
    Vector7d obj_pose;
    bool isBoxConstraints_0 = false;
    bool isBoxConstraints_1 = false;
    Vector3d ub_0;
    Vector3d ub_1;
    Vector3d lb_0;
    Vector3d lb_1;

    HandIK* hand_ik;

    DartDDHand(int fingertype);
    
    DartDDHand(const DartDDHand& dd): DartManipulatorTemplate(dd), 
        radius(dd.radius), fingertype(dd.fingertype),
        finger_length(dd.finger_length), hand_pose(dd.hand_pose){}


    void setHandFrameTransform(const Vector7d& pose);

    void setupCollisionGroup(WorldPtr world);

    void getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips);

    void setConfig(const VectorXd& config, const Vector7d& object_pose);

    bool resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips); 

    void Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, 
                                    std::vector<ContactPoint>* point_contacts);

    bool isIK(const VectorXd& config, const Vector7d& object_pose);

    bool ifIKsolution(const VectorXd& mnp_config, const Vector7d& object_pose){
        return this->isIK(mnp_config, object_pose);
    }

    bool setBoxConstraints(int side, const Vector3d& ub, const Vector3d& lb);

    VectorXd getHandFrameConfig(const VectorXd& config, const Vector7d& object_pose);

};

