#include "cmgmp/manipulators/DartPointManipulator.h"
#include "cmgmp/worlds/DartWorld.h"
#ifndef DART_UTILS
#define DART_UTILS
    #include "cmgmp/dart_utils/dart_utils.h"
#endif
#include "cmgmp/contacts/contact_utils.h"


const double default_shape_density = 1000;  // kg/m^3
const double default_shape_height = 0.1;    // m
const double default_shape_width = 0.03;    // m
const double default_skin_thickness = 1e-3; // m

const double default_start_height = 0.4; // m

const double minimum_start_v = 2.5; // m/s
const double maximum_start_v = 4.0; // m/s
const double default_start_v = 3.5; // m/s

const double minimum_launch_angle = dart::math::toRadian(30.0); // rad
const double maximum_launch_angle = dart::math::toRadian(70.0); // rad
const double default_launch_angle = dart::math::toRadian(45.0); // rad

const double maximum_start_w = 6 * dart::math::constantsd::pi(); // rad/s
const double default_start_w = 3 * dart::math::constantsd::pi(); // rad/s

const double ring_spring_stiffness = 0.5;
const double ring_damping_coefficient = 0.05;
const double default_damping_coefficient = 0.001;

const double default_ground_width = 2;
const double default_wall_thickness = 0.1;
const double default_wall_height = 1;
const double default_spawn_range = 0.9 * default_ground_width / 2;

const double default_restitution = 0.6;

const double default_vertex_stiffness = 1000.0;
const double default_edge_stiffness = 1.0;
const double default_soft_damping = 5.0;

// using namespace dart::dynamics;
// using namespace dart::simulation;
// using namespace dart::collision;
// using namespace dart::common;
// using namespace dart::math;

void setAllColors(const SkeletonPtr& object, const Eigen::Vector3d& color)
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

SkeletonPtr createGround()
{
  SkeletonPtr ground = Skeleton::create("ground");

  BodyNode* bn = ground->createJointAndBodyNodePair<WeldJoint>().second;

  std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(Eigen::Vector3d(
      default_ground_width, default_ground_width, default_wall_thickness));
  auto shapeNode
      = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
          shape);
  shapeNode->getVisualAspect()->setColor(Eigen::Vector3d(1.0, 1.0, 1.0));
  shapeNode->getVisualAspect()->setAlpha(0.5);


  return ground;
}

SkeletonPtr createWall()
{
  SkeletonPtr wall = Skeleton::create("wall");

  BodyNode* bn = wall->createJointAndBodyNodePair<WeldJoint>().second;

  std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(Eigen::Vector3d(
      default_wall_thickness, default_ground_width, default_wall_height));
  auto shapeNode
      = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
          shape);
  shapeNode->getVisualAspect()->setColor(Eigen::Vector3d(0.8, 0.8, 0.8));
  shapeNode->getVisualAspect()->setAlpha(0.5);

  Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
  tf.translation() = Eigen::Vector3d(
      (default_ground_width + default_wall_thickness) / 2.0,
      0.0,
      (default_wall_height - default_wall_thickness) / 2.0);
  bn->getParentJoint()->setTransformFromParentBodyNode(tf);

  shapeNode->getDynamicsAspect()->setRestitutionCoeff(0.2);

  return wall;
}

SkeletonPtr createObject()
{
    SkeletonPtr object = Skeleton::create("rigid_box");

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<FreeJoint>().second;

    // Make the shape based on the requested Shape type
    std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(Eigen::Vector3d(
        default_shape_width, default_shape_width, default_shape_height));

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    // Setup the inertia for the body
    Inertia inertia;
    double mass = default_shape_density * shape->getVolume();
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));
    bn->setInertia(inertia);

    // Set the coefficient of restitution to make the body more bouncy
    shapeNode->getDynamicsAspect()->setRestitutionCoeff(default_restitution);

    setAllColors(object, dart::Color::Red());

    Eigen::Vector6d positions(Eigen::Vector6d::Zero());

    positions[3] = 0.1;
    positions[4] = 0.2;
    positions[5] = default_wall_thickness;
    object->getJoint(0)->setPositions(positions);

    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}

SkeletonPtr createObjectfromMesh(const std::string & filePath)
{
    SkeletonPtr object = Skeleton::create("rigid_ball");

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<FreeJoint>().second;

    // Make the shape based on the requested Shape type
    std::shared_ptr<MeshShape> shape = 
        std::make_shared<MeshShape>(Eigen::Vector3d(3,3,3),MeshShape::loadMesh(filePath));

    shape->setColorMode(MeshShape::ColorMode::MATERIAL_COLOR);

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    // Setup the inertia for the body
    Inertia inertia;
    double mass = default_shape_density * shape->getVolume();
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));
    bn->setInertia(inertia);

    // Set the coefficient of restitution to make the body more bouncy
    // shapeNode->getDynamicsAspect()->setRestitutionCoeff(default_restitution);

    // setAllColors(object, dart::Color::Red());

    Eigen::Vector6d positions(Eigen::Vector6d::Zero());

    positions[3] = 0.1;
    positions[4] = 0.2;
    positions[5] = default_wall_thickness;
    object->getJoint(0)->setPositions(positions);

    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}

SkeletonPtr createCompoundObject(){
    SkeletonPtr object = Skeleton::create("rigid_box");

    // Create the Joint and Body pair
    BodyNode* bn = object->createJointAndBodyNodePair<FreeJoint>().second;

    {
    // Make the shape based on the requested Shape type
    std::shared_ptr<BoxShape> shape = std::make_shared<BoxShape>(Eigen::Vector3d(
        default_shape_width, default_shape_width, default_shape_height));

    auto shapeNode
        = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape);

    // Setup the inertia for the body
    Inertia inertia;
    double mass = default_shape_density * shape->getVolume();
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));
    bn->setInertia(inertia);

    // Set the coefficient of restitution to make the body more bouncy
    shapeNode->getDynamicsAspect()->setRestitutionCoeff(default_restitution);
    }
    // second object
    WeldJoint::Properties properties;
    
    {
      properties.mName = "joint1";
      Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
      tf.translation() = Eigen::Vector3d(0, 0, default_shape_height / 2.0);
      properties.mT_ParentBodyToJoint = tf;
      properties.mT_ChildBodyToJoint = tf.inverse();
    }
    // Create the Joint and Body pair
    BodyNode* bn1 = object->createJointAndBodyNodePair<WeldJoint>(
                         bn, properties,BodyNode::AspectProperties("bodynode1")).second;

    {
      std::shared_ptr<BoxShape> shape1 = std::make_shared<BoxShape>(Eigen::Vector3d(
          default_shape_width, default_shape_width, default_shape_height));
      auto shapeNode1
        = bn1->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(
            shape1);
      
      // Setup the inertia for the body
      Inertia inertia1;
      double mass1 = default_shape_density * shape1->getVolume();
      inertia1.setMass(mass1);
      inertia1.setMoment(shape1->computeInertia(mass1));
      bn1->setInertia(inertia1);

      // Set the coefficient of restitution to make the body more bouncy
      shapeNode1->getDynamicsAspect()->setRestitutionCoeff(default_restitution);
    }
    setBodyNodeColor(bn, dart::Color::Red());
    setBodyNodeColor(bn1, dart::Color::Blue());


    Eigen::Vector6d positions(Eigen::Vector6d::Zero());

    positions[3] = 0.1;
    positions[4] = 0.2;
    positions[5] = default_wall_thickness;
    object->getJoint(0)->setPositions(positions);

    // Add the object to the world
    // object->setName(object->getName() + std::to_string(1));
    return object; 
}

SkeletonPtr createManipulator()
{
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  SkeletonPtr manipulator
      // = loader.parseSkeleton("/home/xianyi/libraries/bullet3/data/kuka_iiwa/model.urdf");
      = loader.parseSkeleton("/home/xianyi/research/CMGMP/data/urdf/KR5/KR5 sixx R650.urdf");
  manipulator->setName("manipulator");

  // Position its base in a reasonable way
  Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
  tf.translation() = Eigen::Vector3d(0.65, 0.0, 0.0);
  manipulator->getJoint(0)->setTransformFromParentBodyNode(tf);

  // Get it into a useful configuration
  // manipulator->getDof(1)->setPosition(dart::math::toRadian(10.0));
  // manipulator->getDof(2)->setPosition(dart::math::toRadian(-10.0));

  std::cout << "Robot Number of DoFs: " << manipulator->getNumDofs() << std::endl;
  for (size_t k = 0; k < manipulator->getNumDofs(); k++){
    manipulator->getDof(k)->setPosition(dart::math::toRadian(10.0));
  }
  std::cout << manipulator->getPositions() << std::endl;
  
  return manipulator;
}

int main(int argc, char* argv[])
{
    DartWorld dw;

    SkeletonPtr ground = createGround();
    SkeletonPtr wall = createWall();

    SkeletonPtr object = createObject();

    // SkeletonPtr robot = createManipulator();
    
    dw.addEnvironmentComponent(ground);
    dw.addEnvironmentComponent(wall);
    
    dw.addObject(object);

    DartPointManipulator* rpt = new DartPointManipulator(2, 0.01);
    
    dw.addRobot(rpt);
    VectorXd rc(12);
    rc << 0,0,0.5,0,0,1,0,0.1,0.5,0,0,1;
    dw.updateRobotConfig(rc);

    std::vector<ContactPoint> pts;
    dw.getObjectContacts(&pts);
    print_contacts(pts);

    // Eigen::VectorXd robot_config(6);
    // robot_config << 0,0,0,0,0,1;
    std::cout << "Robot collision " << dw.isRobotCollide(rc) << std::endl;

    std::vector<Vector7d> object_poses;
    std::vector<VectorXd> mnp_configs;

    Vector7d p1;
    p1 << 0, 0.1, 0.1, 0, 0, 0, 1;
    object_poses.push_back(p1);
    p1[0] = 0.2;
    object_poses.push_back(p1);

    mnp_configs.push_back(rc);
    mnp_configs.push_back(rc);

    dw.setPlaybackTrajectory(object_poses, mnp_configs);
    // dw.setSurfacePoints(surface);

    dw.startWindow(&argc, argv);
    

}