#include <modus/modes/geometry/arrangements.hpp>


int main() {
    int d = 8;
    Eigen::MatrixXd A(d, d);
    A.setIdentity();
    Eigen::VectorXd b(d);

    IncidenceGraphPtr I = initial_arrangement(A, b, 1e-10);

    I.reset();
}