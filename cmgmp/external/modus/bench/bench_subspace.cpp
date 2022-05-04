#include <modus/common/linear_algebra.hpp>
#include <chrono>
#include <iostream>


int main() {
    // Set random seed.
    srand(0);

    // Benchmark null space.
    for (int n = 1; n < 10; n++) {
        for (int m = 2; m < 10; m++) {
            Eigen::MatrixXd A(n, m);
            Eigen::MatrixXd image;
            Eigen::MatrixXd kernel;
            auto start = std::chrono::high_resolution_clock::now();
            for (int k = 0; k < 1000; k++) {
                A.setRandom();
                image_and_kernel_bases(A, image, kernel, 1e-8);
            }
            auto end = std::chrono::high_resolution_clock::now();
            std::cout << n << " " << m << " " <<
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6 / 1e3
            << " ms" << std::endl;
        }
    }
}