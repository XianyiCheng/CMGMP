#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>
// #include <dart/collision/bullet/BulletCollisionDetector.hpp>

#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "../contacts/contact_constraints.h"
#endif

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#ifndef MANIPULATORS_DARTMANIPULATORTEMPLATE
#define MANIPULATORS_DARTMANIPULATORTEMPLATE
    #include "../manipulators/DartManipulatorTemplate.h"
#endif

#ifndef _WORLD_TEMPLATE
#define _WORLD_TEMPLATE
  #include "WorldTemplate.h"
#endif

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::collision;
using namespace dart::common;


#define NOTHING  0
#define PLAY_BACK  1
#define SHOW_OBJECT_ALL  2


class DartWorldWindow : public dart::gui::glut::SimWindow
{
public:
  std::vector<Eigen::Vector6d> object_positions;
  std::vector<Eigen::VectorXd> robot_configs;
  std::vector<ContactPoint> pts;
  SkeletonPtr object;
  DartManipulatorTemplate* robot;
  std::shared_ptr<CollisionGroup> environmentCollisionGroup = 0;
  std::shared_ptr<CollisionGroup> objectCollisionGroup = 0;
  std::shared_ptr<CollisionGroup> robotCollisionGroup = 0;
  std::size_t frameCount = 0;

  int play_mode = NOTHING;
  int mDrawPoints = false;

  DartWorldWindow(
      const WorldPtr& world)
  {
    mBackground[0] = 0.7;
    mBackground[1] = 0.7;
    mBackground[2] = 0.8;
    mBackground[3] = 1;
    setWorld(world);
  }

  void initLights() override
  {
    static float ambient[] = {0.5, 0.5, 0.5, 1.0};
    static float diffuse[] = {0.6, 0.6, 0.6, 1.0};
    static float front_mat_shininess[] = {60.0};
    static float front_mat_specular[] = {0.2, 0.2, 0.2, 1.0};
    static float front_mat_diffuse[] = {0.5, 0.28, 0.38, 1.0};
    static float lmodel_ambient[] = {0.2, 0.2, 0.2, 1.0};
    static float lmodel_twoside[] = {GL_FALSE};

    GLfloat position[] = {1.0, 0.0, 3, 0.0};
    GLfloat position1[] = {-1.0, 0.0, 0, 0.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);

    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, position1);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, front_mat_diffuse);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);
  }

  void keyboard(unsigned char key, int x, int y) override
  {
        if (key == 'a'){
          play_mode = SHOW_OBJECT_ALL;
        }
        SimWindow::keyboard(key, x, y);
  }

  void drawWorld() const override
  {
    // Make sure lighting is turned on and that polygons get filled in
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    SimWindow::drawWorld();
  }

  void draw() override
  {
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    // show object contacts with the environment
    if (mShowMarkers && objectCollisionGroup && environmentCollisionGroup)
    {
      dart::collision::CollisionOption option;
      dart::collision::CollisionResult result;
      objectCollisionGroup->collide(environmentCollisionGroup.get(), option, &result);
      for (const auto& contact : result.getContacts())
      {
        Eigen::Vector3d v = contact.point;
        Eigen::Vector3d f = contact.normal / 50.0;
        glBegin(GL_LINES);
        glVertex3f(v[0], v[1], v[2]);
        glVertex3f(v[0] + f[0], v[1] + f[1], v[2] + f[2]);
        glEnd();
        mRI->setPenColor(Eigen::Vector3d(0.2, 0.8, 0.2));
        mRI->pushMatrix();
        glTranslated(v[0], v[1], v[2]);
        mRI->drawSphere(0.001);
        mRI->popMatrix();
      }
    }
    if (mDrawPoints)
    {
      for (const auto& pt:pts)
      {
        Eigen::Vector3d v = pt.p;
        Eigen::Vector3d f = pt.n/100;
        glBegin(GL_LINES);
        glVertex3f(v[0], v[1], v[2]);
        glVertex3f(v[0] + f[0], v[1] + f[1], v[2] + f[2]);
        glEnd();
        mRI->setPenColor(Eigen::Vector3d(0.2, 0.2, 0.8));
        mRI->pushMatrix();
        glTranslated(v[0], v[1], v[2]);
        mRI->drawSphere(0.001);
        mRI->popMatrix();
      }
    }
    drawWorld();
  }  

  void timeStepping();

};



class DartWorld: public WorldTemplate {
    
    public:
        DartManipulatorTemplate* robot_ptr;
        
        DartWorld();
        DartWorld(const DartWorld& dw): 
              WorldTemplate(dw), world(dw.world), window(dw.window),
              environment(dw.environment), object(dw.object),
              environmentCollisionGroup(dw.environmentCollisionGroup),
              objectCollisionGroup(dw.objectCollisionGroup),
              robotCollisionGroup(dw.robotCollisionGroup){}

        void addEnvironmentComponent(const SkeletonPtr& env);

        void addObject(const SkeletonPtr& the_object);

        void addRobot(DartManipulatorTemplate* the_robot);

        void updateObjectPose(const Vector7d & object_pose);

        void updateRobotConfig(const Eigen::VectorXd & robot_config);

        bool isRobotCollide(const Eigen::VectorXd & robot_config);

        void getObjectContacts(std::vector<ContactPoint>* pts);

        void getObjectContacts(std::vector<ContactPoint>* pts, const Vector7d & object_pose);

        void startWindow(int* pargc, char** argv);

        void setPlaybackTrajectory(const std::vector<Vector7d>& object_traj, const std::vector<VectorXd>& robot_traj);

        void setObjectTrajectory(const std::vector<Vector7d>& object_traj);

        void setSurfacePoints(const std::vector<ContactPoint>& pts);

        ManipulatorTemplate* getRobot(){
          return robot_ptr;
        }


    protected:
        WorldPtr world;
        std::shared_ptr<DartWorldWindow> window;
        std::vector<SkeletonPtr> environment;
        SkeletonPtr object;
        // SkeletonPtr robot;
        std::shared_ptr<CollisionGroup> environmentCollisionGroup = 0; 
        std::shared_ptr<CollisionGroup> objectCollisionGroup = 0;
        std::shared_ptr<CollisionGroup> robotCollisionGroup = 0; 
        

};

