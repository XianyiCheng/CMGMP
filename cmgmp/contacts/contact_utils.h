
#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "contact_kinematics.h"
#endif 
#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "contact_constraints.h"
#endif

// #include <glog/logging.h>
#include <ctime>

#include <iostream>
#include <fstream>  

void surfacepoints_from_file(std::string fileToOpen, std::vector<ContactPoint>* pts);