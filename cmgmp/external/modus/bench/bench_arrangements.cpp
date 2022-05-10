#include <modus/modes/geometry/arrangements.hpp>
#include <chrono>
#include <iostream>
// #include <glog/logging.h>


template <std::size_t alignment>
inline bool is_aligned(void * ptr) noexcept {
    std::size_t max = 1u;
    return std::align(alignment, 1u, ptr, max);
}

int main(int argc, char* argv[]) {
    // google::InitGoogleLogging(argv[0]);

    // Set random seed.
    srand(0);

    // // Benchmark initial arrangements.
    // for (int n = 5; n < 14; n++) {
    //     Eigen::MatrixXd A(n,n);
    //     Eigen::VectorXd b(n);

    //     A.setRandom();
    //     b.setRandom();
    //     auto start = std::chrono::high_resolution_clock::now();
    //     IncidenceGraphPtr I = initial_arrangement(A, b, 1e-8);
    //     auto end = std::chrono::high_resolution_clock::now();
    // }

    // Benchmark increment arrangements.
    for (int j = 0; j < 20; j++)
    {
        int n = 12;
        int d = 11;
        Eigen::MatrixXd A(d,d);
        Eigen::VectorXd b(d);
        A.setRandom();
        b.setRandom();

        IncidenceGraph* I = initial_arrangement(A, b, 1e-8);

        I->_num_arcs_created = 0;

        for (int i = d; i < n; i++) {
            Eigen::VectorXd a(d);
            Eigen::VectorXd b(1);
            a.setRandom();
            b.setRandom();
            increment_arrangement(a, b[0], I, 1e-8);
        }

        std::cout << "num arcs " << I->_num_arcs_created << std::endl;

        delete I;
    }
}