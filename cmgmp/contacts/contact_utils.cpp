#include "contact_utils.h"
#include "../utilities/sample_grasp.h"

void surfacepoints_from_file(std::string fileToOpen, std::vector<ContactPoint>* pts){
    MatrixXd data = load_points_from_csv(fileToOpen);
    int N = data.rows();
    for (int i = 0; i < N; ++i)  {
        ContactPoint pt((data.block(i,0,1,3)).transpose(), -(data.block(i,3,1,3)).transpose(), 0);
        pts->push_back(pt);
    }
}