#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#ifndef MANIPULATORS_MANIPULATORTEMPLATE   
#define MANIPULATORS_MANIPULATORTEMPLATE   
    #include "ManipulatorTemplate.h"
#endif

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::collision;
using namespace dart::common;


class DartManipulatorTemplate: public ManipulatorTemplate
{

public:
    std::vector<SkeletonPtr> bodies;
    std::shared_ptr<CollisionGroup> mCollisionGroup = 0;
    bool ifCheckObjectCollision = false;

    DartManipulatorTemplate(){}
    DartManipulatorTemplate(const DartManipulatorTemplate& dmt): ManipulatorTemplate(dmt), bodies(dmt.bodies), 
                                                                mCollisionGroup(dmt.mCollisionGroup),
                                                                ifCheckObjectCollision(dmt.ifCheckObjectCollision){}
    ~DartManipulatorTemplate(){}

    void addBody(const SkeletonPtr& body){
        bodies.push_back(body);
    }

    virtual void setupCollisionGroup(WorldPtr world) = 0;

    virtual void setConfig(const VectorXd& config, const Vector7d& object_pose) = 0;


};
