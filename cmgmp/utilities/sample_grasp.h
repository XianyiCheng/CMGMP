#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <string>

// void sample_points_box(double hl, double hw, double hh){

// }
using namespace Eigen;
void saveData(std::string fileName, MatrixXd  matrix);
MatrixXd openData(std::string fileToOpen);
MatrixXd load_points_from_csv(std::string fileToOpen);