#include "dart_utils.h"


void setSkeletonColor(const SkeletonPtr& object, const Eigen::Vector3d& color)
{
  // Set the color of all the shapes in the object
  for (std::size_t i = 0; i < object->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = object->getBodyNode(i);
    auto visualShapeNodes = bn->getShapeNodesWith<VisualAspect>();
    for (auto visualShapeNode : visualShapeNodes)
      visualShapeNode->getVisualAspect()->setColor(color);
  }
}

void setBodyNodeColor(BodyNode* bn, const Eigen::Vector3d& color)
{
  // Set the color of all the shapes in the bodynode

  auto visualShapeNodes = bn->getShapeNodesWith<VisualAspect>();
  for (auto visualShapeNode : visualShapeNodes)
    visualShapeNode->getVisualAspect()->setColor(color);

}

SkeletonPtr createFreeBall(const std::string& name, double radius, 
                            const Eigen::Vector3d& color)
{
  SkeletonPtr ball = Skeleton::create(name);

  BodyNode* bn = ball->createJointAndBodyNodePair<FreeJoint>().second;

  std::shared_ptr<EllipsoidShape> shape =std::make_shared<EllipsoidShape>(
      2*radius*Eigen::Vector3d::Ones());
  auto shapeNode
      = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
          shape);
  shapeNode->getVisualAspect()->setColor(color);

  shapeNode->getVisualAspect()->setAlpha(0.85);


  Eigen::Vector6d positions(Eigen::Vector6d::Zero());

  ball->getJoint(0)->setPositions(positions);

  return ball;
}

SkeletonPtr createFreeBox(const std::string& name, const Eigen::Vector3d& dim, 
                            const Eigen::Vector3d& color)
{
    SkeletonPtr object = Skeleton::create(name);

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<FreeJoint>().second;

    // Make the shape based on the requested Shape type
    std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(dim);

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    shapeNode->getVisualAspect()->setColor(color);
    shapeNode->getVisualAspect()->setAlpha(0.85);


    Eigen::Vector6d positions(Eigen::Vector6d::Zero());

    object->getJoint(0)->setPositions(positions);


    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}

SkeletonPtr createFreeCylinder(const std::string& name, double radius, double height, 
                            const Eigen::Vector3d& color)
{
  SkeletonPtr cylinder = Skeleton::create(name);

  BodyNode* bn = cylinder->createJointAndBodyNodePair<FreeJoint>().second;

  std::shared_ptr<CylinderShape> shape =std::make_shared<CylinderShape>(radius, height);
  auto shapeNode
      = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
          shape);
  shapeNode->getVisualAspect()->setColor(color);

  shapeNode->getVisualAspect()->setAlpha(0.85);


  Eigen::Vector6d positions(Eigen::Vector6d::Zero());

  cylinder->getJoint(0)->setPositions(positions);

  return cylinder;
}


SkeletonPtr createFixedBox(const std::string& name, const Eigen::Vector3d& dim, 
                             const Eigen::Vector3d& pos, 
                             const Eigen::Vector3d& color)
{
    SkeletonPtr object = Skeleton::create(name);

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<WeldJoint>().second;

    // Make the shape based on the requested Shape type
    std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(dim);

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    shapeNode->getVisualAspect()->setColor(color);
    shapeNode->getVisualAspect()->setAlpha(0.7);

    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    tf.translation() = pos;
    bn->getParentJoint()->setTransformFromParentBodyNode(tf);

    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}


SkeletonPtr createFreeObjectfromMesh(const std::string & name, const std::string & filePath, 
                                      const Eigen::Vector3d& scale)
{
    SkeletonPtr object = Skeleton::create(name);

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<FreeJoint>().second;

    // Make the shape based on the requested Shape type
    std::shared_ptr<MeshShape> shape = 
        std::make_shared<MeshShape>(scale, MeshShape::loadMesh(filePath));

    shape->setColorMode(MeshShape::ColorMode::MATERIAL_COLOR);

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    // Setup the inertia for the body
    // Inertia inertia;
    // double mass = default_shape_density * shape->getVolume();
    // inertia.setMass(mass);
    // inertia.setMoment(shape->computeInertia(mass));
    // bn->setInertia(inertia);

    // Set the coefficient of restitution to make the body more bouncy
    // shapeNode->getDynamicsAspect()->setRestitutionCoeff(default_restitution);

    // setAllColors(object, dart::Color::Red());

    Eigen::Vector6d positions(Eigen::Vector6d::Zero());
    object->getJoint(0)->setPositions(positions);

    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}



SkeletonPtr createRobot(const std::string & name, const std::string & filePath)
{
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  SkeletonPtr manipulator
      = loader.parseSkeleton(filePath);
  manipulator->setName(name);

  // Position its base in a reasonable way
//   Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
//   tf.translation() = Eigen::Vector3d(0.65, 0.0, 0.0);
//   manipulator->getJoint(0)->setTransformFromParentBodyNode(tf);

//   std::cout << "Robot Number of DoFs: " << manipulator->getNumDofs() << std::endl;
  for (size_t k = 0; k < manipulator->getNumDofs(); k++){
    manipulator->getDof(k)->setPosition(0);
  }
//   std::cout << manipulator->getPositions() << std::endl;
  
  return manipulator;
}

// pose7d: px, py, pz, qx, qy, qz, qw
Vector7d pose6d_to_pose7d(const Eigen::Vector6d & object_pose_6d){
    double q[4];
    so32quat(object_pose_6d.head(3), q);

    Vector7d pose;
    pose.head(3) = object_pose_6d.tail(3);
    pose[3] = q[1];
    pose[4] = q[2];
    pose[5] = q[3];
    pose[6] = q[0];

    return pose;

}

Eigen::Vector6d pose7d_to_pose6d(const Vector7d & pose_7d){
    Eigen::Quaterniond q(pose_7d[6], pose_7d[3], pose_7d[4], pose_7d[5]);
    Eigen::AngleAxisd aa(q);
    Eigen::Vector3d aaxis = aa.axis()*aa.angle();
    Eigen::Vector6d pose6d;
    pose6d.head(3) = aaxis;
    pose6d.tail(3) = pose_7d.head(3);

    return pose6d;
    
}