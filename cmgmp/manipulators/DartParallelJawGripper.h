#ifndef MANIPULATORS_DARTMANIPULATORTEMPLATE
#define MANIPULATORS_DARTMANIPULATORTEMPLATE
    #include "DartManipulatorTemplate.h"
#endif

#ifndef MANIPULATORS_MANIPULATORTEMPLATE   
#define MANIPULATORS_MANIPULATORTEMPLATE   
    #include "ManipulatorTemplate.h"
#endif

#ifndef DART_UTILS
#define DART_UTILS
    #include "../dart_utils/dart_utils.h"
#endif

class DartParallelJawGripper: public virtual DartManipulatorTemplate {

public:
    // config: [p1, n1, p2, n2, n_finger] or [p1, n1, n_gripper, n_finger]

    double fingertip_radius = 0;
    double max_open_dist = 0;
    double finger_length = 0;

    DartParallelJawGripper(double max_d, double radius, double fl = 0);
    DartParallelJawGripper(const DartParallelJawGripper& dpm): DartManipulatorTemplate(dpm), fingertip_radius(dpm.fingertip_radius), max_open_dist(dpm.max_open_dist){}

    void setupCollisionGroup(WorldPtr world);

    void getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips);

    void setConfig(const VectorXd& config, const Vector7d& object_pose);

    bool resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface, 
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips); 

    void Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, 
                                    std::vector<ContactPoint>* point_contacts);

    
    bool ifIK(const VectorXd& config);

    bool ifIK(const std::vector<ContactPoint>& fingertips);


    

};