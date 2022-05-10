#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::collision;
using namespace dart::common;

void setSkeletonColor(const SkeletonPtr& object, const Eigen::Vector3d& color);

void setBodyNodeColor(BodyNode* bn, const Eigen::Vector3d& color);

SkeletonPtr createFreeBall(const std::string& name, double radius, 
                            const Eigen::Vector3d& color = Eigen::Vector3d(0.7,0.3,0.3));

SkeletonPtr createFreeBox(const std::string& name, const Eigen::Vector3d& dim, 
                            const Eigen::Vector3d& color = Eigen::Vector3d(0.7,0.3,0.3));

SkeletonPtr createFreeCylinder(const std::string& name, double radius, double height, 
                            const Eigen::Vector3d& color = Eigen::Vector3d(0.7,0.3,0.3));

SkeletonPtr createFixedBox(const std::string& name, const Eigen::Vector3d& dim, 
                            const Eigen::Vector3d& pos, const Eigen::Vector3d& color = Eigen::Vector3d(0.4,0.4,0.4));

SkeletonPtr createFreeObjectfromMesh(const std::string & name, const std::string & filePath, 
                                        const Eigen::Vector3d& scale = Eigen::Vector3d(1,1,1));

SkeletonPtr createRobot(const std::string & name, const std::string & filePath);

Vector7d pose6d_to_pose7d(const Eigen::Vector6d & object_pose_6d);

Eigen::Vector6d pose7d_to_pose6d(const Vector7d & pose_7d);
