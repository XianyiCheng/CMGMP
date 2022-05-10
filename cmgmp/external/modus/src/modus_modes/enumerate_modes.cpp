#include <modus/common/linear_algebra.hpp>
#include <modus/modes/enumerate.hpp>
#include <modus/modes/geometry/arrangements.hpp>
#include <modus/modes/geometry/interior_point.hpp>
#include <modus/modes/geometry/convex_hull.hpp>
#include <modus/modes/geometry/halfspace_intersection.hpp>
#include <modus/modes/preprocess/matroid_preprocessor.hpp>
#include <modus/common/assert.hpp>
// #include <modus/common/logging.hpp>
// #include <glog/logging.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>

static int DEBUG=0;
static int PROFILE=0;

using namespace modus;

IncidenceGraph* enumerate_cs_modes(const Eigen::MatrixXd& N, 
                                   const Eigen::VectorXd& b, double eps,
                                   ModeEnumerationOptions* options) {
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
    if (options && options->interior_point.size() > 0) {
        r = options->interior_point;
    } else {
        r = interior_point(N, b, eps);
    }
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  int pt: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

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

    // Get affine subspace.
    Eigen::MatrixXd affine = orth(N_s.transpose(), eps);
    // Project hyperplanes and interior point into affine subspace.
    N_s = N_s * affine;
    N_c = N_c * affine;
    N_proj = N_proj * affine;
    r = affine.transpose() * r;

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
    // if (DEBUG) {
    //     DLOG(INFO) << "\n\t" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << "dual points:\n" << dual;
    // }

    // Remove duplicate points.
    Eigen::VectorXi idx_unique = UniqueRows(dual_aff.transpose(), eps);
    // std::cout << "idx unique " << idx_unique.transpose() << std::endl;
    dual_aff = get_cols(dual_aff, idx_unique);
    N_s = GetRows(N_s, idx_unique);
    b_s = GetRows(b_s, idx_unique);
    n_s = N_s.rows();
    // if (DEBUG) {
    //     DLOG(INFO) << "N_s\n" << std::fixed << std::setprecision(6) << std::setfill(' ') 
    //     << N_s << std::endl;
    // }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "partproj: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    // TODO Filter duplicate points.
    // DLOG(INFO) << std::fixed << std::setfill(' ') << 
    //  "\ndual_aff\n" << dual_aff << std::endl;
    // // TODO Handle additional degenerate cases.

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
    // Update interior points.
    graph->set_hyperplanes(N_s, b_s);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->update_interior_point(eps);
        }
    }

    // Lift interior points and update sign vectors.
    graph->set_hyperplanes(N, b);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->interior_point = kernel * affine * u->interior_point;
        }
    }
    graph->update_sign_vectors(eps);

    // Mark all proper faces as feasible.
    for (int k = 0; k < d_eff + 2; k++) {
        for (Node* u : graph->rank(k)) {
            u->feasible = true;
        }
    }

    // std::cout << "HEREHERHERE"<< std::endl;
    // for (Node* u : graph->_nodes) {
    //     std::cout << *u << std::endl;
    // }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
            << " ms" << std::endl << std::endl;
    }

    return graph;
}

IncidenceGraph* enumerate_ss_modes(const Eigen::MatrixXd& N,
                                   const Eigen::VectorXd& d,
                                   const Eigen::MatrixXd& T, 
                                   std::string cs_mode, double eps,
                                   ModeEnumerationOptions* options) {
    // Algorithm steps:
    // 1. Partition the hyperplanes into N_c and [N_s|T_c].
    // 2. Project into subspaces of N_c and [N_s|T_c].
    //      2.a. Null space of N_c,
    //      2.b. Row space of [N_s|T_c].
    // 3. Pre-process hyperplanes.
    //      3.a. Remove non-unique rows.
    //      3.b. Re-order for linear independence (if using 4.b.).
    //      3.c. Remove zero rows.
    // 4. Initialize the hyperplane arrangement using one of the following:
    //      4.a. Convex hull of N_s,
    //      4.b. Reduced combination of N_s and T_c.
    // 5. Increment the hyperplane arrangement with the remaining hyperplanes.
    // 6. Update positions and sign vectors.

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        std::cout << "ss modes: " << cs_mode << std::endl;
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 1. Partition the hyperplanes into N_c and [N_s|T_c]                    *
     **************************************************************************/
    // Create index vectors.
    Eigen::VectorXi index_Nc, index_Ns, index_Tc;
    arg_where(cs_mode, '0', index_Nc);
    arg_where(cs_mode, '-', index_Ns);
    // Handle alternate format with 'c' and 's'.
    if (index_Nc.size() == 0 && index_Ns.size() == 0) {
        arg_where(cs_mode, 'c', index_Nc);
        arg_where(cs_mode, 's', index_Ns);
    }
    int n_sliding_dir = T.rows() / N.rows();
    assert(T.rows() % N.rows() == 0);
    int n_Tc = 0;
    for (int i = 0; i < T.rows(); i++) {
        int m = i / n_sliding_dir;
        if (cs_mode[m] == '0' || cs_mode[m] == 'c') {
            n_Tc += 1;
        }
    }
    index_Tc.resize(n_Tc);
    int j = 0;
    for (int i = 0; i < T.rows(); i++) {
        int m = i / n_sliding_dir;
        if (cs_mode[m] == '0' || cs_mode[m] == 'c') {
            index_Tc[j++] = i;
        }
    }
    // Partition hyperplanes.
    Eigen::MatrixXd Nc, Ns, Tc, NsTc, N_proj, Tc_raw;
    Eigen::VectorXd dc, ds, b;
    GetRows(N, index_Nc, Nc);
    GetRows(d, index_Nc, dc);
    GetRows(N, index_Ns, Ns);
    GetRows(d, index_Ns, ds);
    GetRows(T, index_Tc, Tc);
    GetRows(T, index_Tc, Tc_raw);
    int dim = N.cols();
    int n_ns = Ns.rows();
    int n_tc = Tc.rows();
    NsTc.resize(n_ns + n_tc, dim);
    NsTc.topRows(n_ns) = Ns;
    NsTc.bottomRows(n_tc) = Tc;
    N_proj = N;
    b.resize(NsTc.rows());
    b.setZero();
    b.block(0,0,n_ns,1) = GetRows(d, index_Ns);
    // if (DEBUG) {
    //     DLOG(INFO) << "\nPartition hyperplanes:\n"
    //                << "\tNc: " << Nc.rows() << "x" << Nc.cols() << "\n"
    //                << "\tNs: " << Ns.rows() << "x" << Ns.cols() << "\n"
    //                << "\tTc: " << Tc.rows() << "x" << Tc.cols() << "\n"
    //                << "\tNsTc: " << NsTc.rows() << "x" << NsTc.cols();
    // }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "    part: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
            << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 2. Project into subspaces of N_c and [N_s|T_c].                        *
     * 2.a. Null space of N_c,                                                *
     * 2.b. Row space of [N_s|T_c].                                           *
     **************************************************************************/
    // Project into null space of N_c.
    Eigen::MatrixXd kernel;
    kernel_basis(Nc, kernel, eps);
    if (kernel.cols() > 0) {
        Ns = Ns * kernel;
        Tc = Tc * kernel;
        NsTc = NsTc * kernel;
        N_proj = N_proj * kernel;
    }
    // Project into row space of [N_s|T_c].
    Eigen::MatrixXd rowspace;
    orth(NsTc.transpose(), eps, rowspace);
    Ns = Ns * rowspace;
    Tc = Tc * rowspace;
    NsTc = NsTc * rowspace;
    N_proj = N_proj * rowspace;
    // if (DEBUG) {
    //     DLOG(INFO) << "\n"
    //     << "Project hyperplanes:\n"
    //     << "\t  Nc null: " << dim << " -> " << kernel.cols() << "\n"
    //     << "\tNsTc orth: " << kernel.cols() << " -> " << rowspace.cols() << "\n"
    //     << PrettyMatrix(Ns, 
    //        "\t       Ns:") << "\n"
    //     << PrettyMatrix(Tc, 
    //        "\t       Tc:") << "\n";
    // }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "    proj: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
            << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 3. Pre-process hyperplanes.                                            *
     *      3.a. Remove non-unique rows.                                      *
     *      3.b. Re-order for linear independence (if using 4.b.).            *
     *      3.c. Remove zero rows.                                            *
     **************************************************************************/
    // Annotate each hyperplane with # of sides.
    Eigen::VectorXi sides(NsTc.rows());
    sides.setOnes();
    sides.topRows(n_ns) *= 2;
    sides.bottomRows(NsTc.rows() - n_ns) *= 3;

    // Normalize hyperplanes.
    for (int i = 0; i < NsTc.rows(); i++) {
        double norm = NsTc.row(i).norm();
        NsTc.row(i) /= norm;
        b[i] /= norm;
    }

    // Remove duplicate hyperplanes, favoring normals (sides = 2) over tangents
    // (sides = 3).
    // Eigen::VectorXi idx = unique_rows(NsTc, eps);
    Eigen::VectorXi idx = UniqueRowsReoriented(NsTc, eps);
    NsTc = GetRows(NsTc, idx);
    sides = GetRows(sides, idx);
    b = GetRows(b, idx);

    // Reorder hyperplanes so that normals (sides = 2) are first.
    idx = LexicographicArgSortRows(sides);
    NsTc = GetRows(NsTc, idx);
    sides = GetRows(sides, idx);
    b = GetRows(b, idx);

    // TODO Reorder hyperplanes to be linearly independent in top d rows.

    // Remove 0 hyperplanes and output.
    idx = NonzeroRows(NsTc, eps);
    NsTc = GetRows(NsTc, idx);
    sides = GetRows(sides, idx);
    b = GetRows(b, idx);

    // Apply matroid preprocessing.
    MatroidPreprocessor matroid_pre;
    matroid_pre.SetVerbose(false);
    matroid_pre.Preprocess(NsTc, eps, 1e2);
    NsTc = matroid_pre.GetMinimallyProjectivelyEquivalentMatrix();
    idx = matroid_pre.GetIndices();
    sides = GetRows(sides, idx);
    b = GetRows(b, idx);

    // Compute dimension of matrices.
    int d_eff = NsTc.cols();
    int d_ns = 0;
    if (Ns.rows() > 0) {
        Eigen::MatrixXd R;
        orth(Ns.transpose(), eps, R);
        d_ns = R.cols();
    }
    n_ns = 0;
    for (int i = 0; i < sides.size(); i++) {
        n_ns += sides[i] == 2;
    }

    // if (DEBUG) {
    //     DLOG(INFO) << "\n"
    //     << "Preprocess hyperplanes:" << "\n"
    //     << "\tNsTc dim: " << NsTc.rows() << "x" << NsTc.cols() << "\n"
    //     << "\t  Ns dim: " << n_ns << "x" << d_ns << "\n"
    //     << PrettyMatrix(NsTc, 
    //        "\t    NsTc:") << "\n"
    //     << PrettyMatrix(b,    
    //        "\t       b:") << "\n";
    // }

    /**************************************************************************
     * 4. Initialize the hyperplane arrangement using one of the following:   *
     *      4.a. Convex hull of N_s,                                          *
     *      4.b. Reduced combination of N_s and T_c.                          *
     **************************************************************************/
    // Create intialization variables.
    IncidenceGraph* graph = nullptr;
    int start_incr = 0;
    bool init_convex = false;

    // If # separating > dim(Ns) and dim(Ns) == dim(NsTc), then initialize with
    // convex hull.
    if (n_ns > d_ns && d_ns == d_eff && d_ns != 1) {
        // Get preprocessed Ns.
        idx = ArgWhereEqual(sides, 2, eps);
        Ns = GetRows(NsTc, idx);
        Eigen::VectorXd bs = GetRows(b, idx);

        // Compute halfspace intersection.
        init_convex = true;
        if (options && options->interior_point.size() > 0) {
            Eigen::MatrixXd L = kernel * rowspace;
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(L);
            Eigen::VectorXd int_pt = qr.solve(options->interior_point);
            graph = halfspace_intersection(Ns, bs, eps, &int_pt);
        } else {
            graph = halfspace_intersection(Ns, bs, eps);
        }

        start_incr = n_ns;
    }
    // Otherwise, intialize with enumeration.
    else {
        // Preprocess hyperplanes, i.e. remove duplicates and re-order.
        partial_preprocess_hyperplanes(NsTc, b, sides, eps);

        // Create initial arrangement.
        graph = partial_initial_arrangement(NsTc, b, sides, eps);
        // graph = initial_arrangement(NsTc.topRows(d_eff), b.topRows(d_eff), eps);
        start_incr = d_eff;
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        if (init_convex) {
        std::cout << "  convex: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        } else {
        std::cout << " partial: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
            << " ms" << std::endl;
        }
        start = std::chrono::high_resolution_clock::now();
    }

    // Exit early if there are no active contacts.
    if (Nc.rows() == 0) {
        if (PROFILE) {
            auto end = std::chrono::high_resolution_clock::now();
            std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
                << " ms" << std::endl;
            start = std::chrono::high_resolution_clock::now();
        }
        return graph;
    }
    
    /**************************************************************************
     * 5. Increment the hyperplane arrangement with the remaining planes.     *
     **************************************************************************/
    for (int i = start_incr; i < NsTc.rows(); i++) {
        increment_arrangement(NsTc.row(i), 0, graph, eps);
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "    incr: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
            << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 6. Update positions and sign vectors.                                  *
     **************************************************************************/
    Eigen::MatrixXd NTc(N_proj.rows() + Tc.rows(), d_eff);
    NTc.topRows(N_proj.rows()) = N_proj;
    NTc.bottomRows(Tc.rows()) = Tc;
    graph->A = NTc;
    graph->b = Eigen::VectorXd::Zero(NTc.rows());
    graph->update_sign_vectors(eps);

    // Lift interior points.
    NTc.resize(N.rows() + Tc_raw.rows(), N.cols());
    NTc.topRows(N.rows()) = N;
    NTc.bottomRows(Tc_raw.rows()) = Tc_raw;
    graph->A = NTc;
    graph->b = Eigen::VectorXd::Zero(NTc.rows());
    for (int k = 0; k <= d_eff; k++) {
        for (Node* u : graph->rank(k)) {
            u->interior_point = kernel * rowspace * u->interior_point;
        }
    }
    graph->update_sign_vectors(eps);

    // Mark nodes as feasible (a valid ss-mode) or not.
    int n_contacts = N.rows();
    for (int k = 0; k <= d_eff; k++) {
        for (Node* u : graph->rank(k)) {
            u->feasible = u->sign_vector.substr(0, n_contacts) == cs_mode;
        }
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  update: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
            << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
            << " ms" << std::endl << std::endl;
    }

    return graph;
}

std::vector<IncidenceGraph*> 
enumerate_ss_modes(IncidenceGraph* cs_graph, const Eigen::MatrixXd& N, 
                   const Eigen::VectorXd& d, const Eigen::MatrixXd& T, double eps,
                   ModeEnumerationOptions* options) {
    // 
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (true) {
        std::cout << "ss modes: " << N.rows() << "+" << T.rows() << "x" << N.cols() << std::endl;
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    std::vector<std::string> cs_modes = cs_graph->get_proper_sign_vectors();
    std::vector<IncidenceGraph*> ss_graphs;
    for (int i = 0; i < cs_modes.size(); i++) {
        // Skip all-separating mode.
        if (cs_modes[i].find_first_not_of('-') == std::string::npos) {
            continue;
        }
        ss_graphs.push_back(enumerate_ss_modes(N, d, T, cs_modes[i], eps, options));
    }

    if (true) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
            << " ms" << std::endl << std::endl;
    }

    return ss_graphs;
}


IncidenceGraph* modus::EnumerateCSModes(const Eigen::MatrixXd& N, 
                                        double eps, 
                                        ModeEnumerationOptions* options)
{
    // We always use 0 contact distance.
    Eigen::VectorXd d = Eigen::VectorXd::Zero(N.rows());

    // Wrap call in try-catch to enable saving the inputs on failure.
    try {
        return enumerate_cs_modes(N, d, eps, options);
    } catch (ppk::assert::AssertionException& e) {
        using namespace modus;

        // std::cout << "Caught exception" << std::endl;

        // Log error.
        // LOG(ERROR) << e.what();

        // Log inputs.
        // WriteInputToLog(CreateInput(N, d, eps));

        // Escalate exception.
        // throw e;
    } catch (...) {
        std::cout << "WTF"<< std::endl;
        // Log inputs.
        // WriteInputToLog(CreateInput(N, d, eps));
    }

    return nullptr;
}

IncidenceGraph* modus::EnumerateSSModes(const Eigen::MatrixXd& N, 
                                        const Eigen::MatrixXd& T,
                                        const std::string& cs_mode, 
                                        double eps, 
                                        ModeEnumerationOptions* options)
{
    // We always use 0 contact distance.
    Eigen::VectorXd d = Eigen::VectorXd::Zero(N.rows());

    // Wrap call in try-catch to enable saving the inputs on failure.
    try {
        return enumerate_ss_modes(N, d, T, cs_mode, eps, options);
    } catch (ppk::assert::AssertionException& e) {
        using namespace modus;

        // std::cout << "Caught exception" << std::endl;

        // Log error.
        // LOG(ERROR) << e.what();

        // Log inputs.
        // WriteInputToLog(CreateInput(N, d, T, cs_mode, eps));

        // Escalate exception.
        // throw e;
    } catch (...) {
        std::cout << "WTF"<< std::endl;
        // Log inputs.
        // WriteInputToLog(CreateInput(N, d, T, cs_mode, eps));
    }

    return nullptr;
}
