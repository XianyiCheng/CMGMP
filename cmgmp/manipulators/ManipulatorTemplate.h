
#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "../contacts/contact_constraints.h"
#endif

class ManipulatorTemplate
{
public:
    int n_pts = 0;

// public:
    ManipulatorTemplate(){};
    ~ManipulatorTemplate(){};

    ManipulatorTemplate(const ManipulatorTemplate& mt): n_pts(mt.n_pts){}

    virtual void setConfig(const VectorXd& config, const Vector7d& object_pose) {};

    virtual bool resampleFingers(int n_on, const VectorXd& config, const Vector7d& object_pose, const std::vector<ContactPoint>& object_surface,
    VectorXd& new_config, std::vector<ContactPoint>* remain_fingertips) =0; 

    virtual void getFingertipsOnObject(const VectorXd& config, const Vector7d& object_pose, std::vector<ContactPoint>* fingertips) = 0;

    // contact model
    virtual void Fingertips2PointContacts(const std::vector<ContactPoint>& fingertips, std::vector<ContactPoint>* point_contacts) =0;

    virtual bool isSameConfig(const VectorXd& config0, const VectorXd& config1){
        if (config0.size()!=config1.size()){
            return false;
        } else {
            if ((config0 - config1).norm() < 1e-3){
                return true;
            } else {
                return false;
            }
        }
    }

    int getNumberofFingertips(){
        return this->n_pts;
    };

    virtual bool ifIKsolution(const VectorXd& mnp_config, const Vector7d& object_pose){
        return true;
    }
    
};