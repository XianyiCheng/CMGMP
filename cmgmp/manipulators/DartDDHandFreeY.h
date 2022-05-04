#include "DartDDHand.h"

class DartDDHandFreeY: public DartDDHand {

public:

    // config: p1, n1, p2, n2, y
    double max_dy = 3E-3;

    DartDDHandFreeY(int fingertype) : DartDDHand(fingertype) {
    }

    void setConfig(const VectorXd& config, const Vector7d& object_pose);

    bool resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips);

    VectorXd getHandFrameConfig(const VectorXd& config, const Vector7d& object_pose); // return [p1, p2, y_world]
};