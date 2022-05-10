#include <modus/modes/geometry/halfspace_intersection.hpp>
#include <modus/modes/geometry/interior_point.hpp>
#include <modus/modes/geometry/convex_hull.hpp>
#include <modus/common/linear_algebra.hpp>
// #include <glog/logging.h>
#include <iomanip>
#include <chrono>
#include <iostream>


static int DEBUG=0;
static int PROFILE=0;

using namespace modus;

IncidenceGraph* halfspace_intersection(const Eigen::MatrixXd& N,
                                       const Eigen::VectorXd& b,
                                       double eps,
                                       const Eigen::VectorXd* int_pt) {
    // 1. Find a strictly interior point (if it exists).
    // 2. Partition hyperplanes.
    // 3. Project into null space, affine space, and dual space. 
    // 4. Build convex hull and face lattice.
    // 5. Update interior points.
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        std::cout << "cs modes: " << N.rows() << "x" << N.cols() << std::endl;
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 1. Find a strictly interior point.                                     *
     **************************************************************************/
    // Get dimensions.
    int n = N.rows();
    int d = N.cols();
    // Find a strictly interior point.
    Eigen::VectorXd r;
    if (!int_pt) {
        r = interior_point(N, b, eps);
    } else {
        r = *int_pt;
    }
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  int pt: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    // if (DEBUG) {
    //     DLOG(INFO) << "\t\nint pt\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << r << std::endl;
    // }

    /**************************************************************************
     * 2. Partition hyperplanes.                                              *
     **************************************************************************/
    // Get partition variables.
    Eigen::MatrixXd N_c, N_s, N_proj;
    Eigen::VectorXd b_c, b_s;
    Eigen::VectorXi idx_c, idx_s;
    int n_s, n_c;
    // Find hyperplanes which are always in contact.
    Eigen::VectorXd res = N*r - b;
    arg_where_zero(res, eps, idx_c);
    arg_where_nonzero(res, eps, idx_s);
    n_c = idx_c.size();
    n_s = idx_s.size();
    // Partition hyperplanes.
    GetRows(N, idx_c, N_c);
    GetRows(N, idx_s, N_s);
    GetRows(b, idx_c, b_c);
    GetRows(b, idx_s, b_s);

    /**************************************************************************
     * 3. Project into null space, affine space, and dual affine space.       *
     **************************************************************************/
    // Get nullspace. 
    Eigen::MatrixXd kernel;
    kernel_basis(N_c, kernel, eps);
    // Project hyperplanes and interior point into contacting nullspace.
    N_s = N_s * kernel;
    N_c = N_c * kernel;
    N_proj = N * kernel;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.setThreshold(eps);
    qr.compute(kernel);
    Eigen::VectorXd tmp = qr.solve(r); // Project interior point.
    r = tmp;
    // if (DEBUG) {
    //     DLOG(INFO) << "\t\nkernel\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << kernel << std::endl;
    //     DLOG(INFO) << "\t\nN_proj\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << N_proj << std::endl;
    // }

    // Get affine subspace.
    Eigen::MatrixXd affine = orth(N_s.transpose(), eps);
    // Project hyperplanes and interior point into affine subspace.
    N_s = N_s * affine;
    N_c = N_c * affine;
    N_proj = N_proj * affine;
    r = affine.transpose() * r;
    // if (DEBUG) {
    //     DLOG(INFO) << "\t\notho\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << affine << std::endl;
    //     DLOG(INFO) << "\t\nN_proj\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << N_proj << std::endl;
    // }

    // Get dual points.
    Eigen::VectorXd b_off = b_s - N_s * r;
    Eigen::MatrixXd dual = N_s;
    for (int i = 0; i < N_s.rows(); i++) {
        dual.row(i) /= b_off[i];
    }

    // Project dual points into affine subspace.
    Eigen::MatrixXd dual_aff;
    project_affine(dual.transpose(), eps, dual_aff);
    int d_eff = dual_aff.rows();
    // DLOG(INFO) << "\n\t" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    // << "dual points:\n" << dual;

    // Remove duplicate points.
    Eigen::VectorXi idx_unique = UniqueRows(dual_aff.transpose(), eps);
    // std::cout << "idx unique " << idx_unique.transpose() << std::endl;
    dual_aff = get_cols(dual_aff, idx_unique);
    N_s = GetRows(N_s, idx_unique);
    b_s = GetRows(b_s, idx_unique);
    n_s = N_s.rows();
    // DLOG(INFO) << "N_s\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    // << N_s << std::endl;

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "partproj: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    // TODO Filter duplicate points.
    // DLOG(INFO) << std::fixed << std::setfill(' ') << 
    //  "\ndual_aff\n" << dual_aff << std::endl;
    // TODO Handle additional degenerate cases.

    /**************************************************************************
     * 4. Build convex hull and face lattice.                                 *
     **************************************************************************/
    // Convex hull.
    Eigen::MatrixXi M = convex_hull(dual_aff, eps);
    // std::cout << "M\n" << std::fixed << M << std::endl;
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "convhull: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    // Build face lattice.
    IncidenceGraph* graph = build_graph_cone(M.transpose(), n_s, d_eff, eps);
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << " lattice: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 5. Update interior point and sign vectors.                             *
     **************************************************************************/
    // Compute interior points.
    graph->set_hyperplanes(N_s, b_s);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->update_interior_point(eps);
            // if (k <= 1) {
            //     std::cout << *u << std::endl;
            // }
            // u->update_sign_vector(eps);
            // assert(u->sign_vector.find('+') == std::string::npos);
        }
    }
    // graph->update_sign_vectors(eps);
    // std::cout << "initial sign vectors" << std::endl;
    // for (int k = 0; k <= d_eff + 1; k++) {
    //     for (Node* u : graph->rank(k)) {
    //         std::cout << u->_key << std::endl;
    //         std::cout << u->sign_vector << std::endl;
    //         std::cout << N_s * u->interior_point << std::endl;
    //         std::cout << std::endl;
    //     }
    // }
    // for (std::string s : graph->get_sign_vectors()) {
    //     std::cout << s << std::endl;
    //     assert(s.find('+') == std::string::npos);
    // }
    
    // std::cout << "update interior point" << std::endl;
    // for (std::string s : graph->get_sign_vectors()) {
    //     std::cout << s << std::endl;
    //     assert(s.find('+') == std::string::npos);
    // }

    // Lift interior points and update sign vectors.
    // std::cout << "lift interior point" << std::endl;
    graph->set_hyperplanes(N, b);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->interior_point = kernel * affine * u->interior_point;
        }
    }
    graph->update_sign_vectors(eps);

    // for (int k = 0; k <= d_eff + 1; k++) {
    //     for (Node* u : graph->rank(k)) {
    //         std::cout << *u << std::endl;
    //     }
    // }
    
    // Update internal keys.
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->_key = u->sign_vector;
        }
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
            << " ms" << std::endl << std::endl;
    }

    return graph;
}