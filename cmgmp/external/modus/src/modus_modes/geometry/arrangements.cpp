#include <modus/modes/geometry/arrangements.hpp>
#include <modus/modes/geometry/interior_point.hpp>
#include <modus/common/linear_algebra.hpp>
#include <modus/modes/geometry/convex_hull.hpp>
#include <modus/common/assert.hpp>
// #include <glog/logging.h>
#include <chrono>
#include <iostream>
#include <list>
#include <iomanip>


static int DEBUG=0;
static int PROFILE=0;

using namespace modus;

// TODO Replace me
int get_position(Node* f, const Eigen::VectorXd& a, double b, double eps) {
    return get_position(a.dot(f->interior_point) - b, eps);
}

void get_vector(Node* e, Eigen::VectorXd& v_e) {
    MODUS_ASSERT(e->rank == 1);
    IncidenceGraph* I = e->_graph;
    v_e.resize(I->A.cols());
    v_e.setZero();
    for (Node* u : e->subfaces) {
        v_e += u->interior_point;
    }
    if (e->subfaces.size() == 2) {
        v_e /= 2.0;
    } else if (e->subfaces.size() == 1) {
        v_e = e->interior_point - v_e;
    } else {
        MODUS_ASSERT(false);
    }
}

int get_color_vertex(Node* v, const Eigen::VectorXd& a, double b, double eps) {
    MODUS_ASSERT(v->rank == 0);
    int s = get_position(v, a, b, eps);
    if (s == 0) {
        return COLOR_AH_CRIMSON;
    } else {
        return COLOR_AH_WHITE;
    }
}

int get_color_edge(Node* e, const Eigen::VectorXd& a, double b, double eps) {
    MODUS_ASSERT(e->rank == 1);
    IncidenceGraph* I = e->_graph;
    if (e->subfaces.size() == 2) {
        auto it = e->subfaces.begin();
        Node* v0 = *it++;
        Node* v1 = *it;
        MODUS_ASSERT(v0 != v1);
        int s0 = get_position(v0, a, b, eps);
        int s1 = get_position(v1, a, b, eps);
        if (DEBUG) {
            // printf("edge pos: %d %d\n", s0, s1);
        }
        if (s0 * s1 == 1) {
            return COLOR_AH_WHITE;
        } else if (s0 == 0 && s1 == 0) {
            return COLOR_AH_CRIMSON;
        } else if (s0 == 0 || s1 == 0) {
            return COLOR_AH_PINK;
        } else if (s0 + s1 == 0) {
            return COLOR_AH_RED;
        } else {
            MODUS_ASSERT(false);
        }
    } else if (e->subfaces.size() == 1) {
        Node* v0 = *e->subfaces.begin();
        int s0 = get_position(v0, a, b, eps);
        Eigen::VectorXd v_e(e->interior_point.size());
        get_vector(e, v_e);
        int s_e = get_position(a.dot(v_e), eps);
        if (DEBUG) {
            // printf("vert pos edge pos: %d %d\n", s0, s_e);
        }
        if (s0 == 0 && s_e == 0) {
            return COLOR_AH_CRIMSON;
        } else if (s0 == 0 && s_e != 0) {
            return COLOR_AH_PINK;
        } else if  (s0 * s_e == 1) {
            return COLOR_AH_WHITE;
        } else if (s0 * s_e == -1) {
            return COLOR_AH_RED;
        } else if (s0 != 0 and s_e == 0) {
            return COLOR_AH_WHITE;
        } else {
            MODUS_ASSERT(false);
        }
    } else {
        MODUS_ASSERT(false);
    }
    return COLOR_AH_WHITE;
}

int get_side_vertex(Node* v, const Eigen::VectorXd& a, double b, double eps) {
    MODUS_ASSERT(v->rank == 0);
    int s = get_position(v, a, b, eps);
    if (s == 0) {
        return SIDE_AH_CRIMSON;
    } else if (s == -1) {
        return SIDE_AH_BELOW;
    } else if (s == 1) {
        return SIDE_AH_ABOVE;
    } else {
        MODUS_ASSERT(false);
    }
    return 0;
}

int get_side_edge(Node* e, const Eigen::VectorXd& a, double b, double eps) {
    MODUS_ASSERT(e->rank == 1);
    IncidenceGraph* I = e->_graph;
    int side = SIDE_AH_UNSET;
    if (e->subfaces.size() == 2) {
        auto it = e->subfaces.begin();
        Node* v0 = *it++;
        Node* v1 = *it;
        MODUS_ASSERT(v0 != v1);
        int s0 = get_position(v0, a, b, eps);
        int s1 = get_position(v1, a, b, eps);
        if (s0 == 0 && s1 == 0) {
            side = SIDE_AH_CRIMSON;
        }
        int s = s0 + s1;
        if (s < 0) {
            side = SIDE_AH_BELOW;
        } else if (s > 0) {
            side = SIDE_AH_ABOVE;
        } else {
            side = SIDE_AH_RED;
        }
    } else if (e->subfaces.size() == 1) {
        Node* v0 = *e->subfaces.begin();
        int s0 = get_position(v0, a, b, eps);
        Eigen::VectorXd v_e(e->interior_point.size());
        get_vector(e, v_e);
        int s_e = get_position(a.dot(v_e), eps);
        if (s0 == 0 && s_e == 0) {
            side = SIDE_AH_CRIMSON;
        } else if (s0 == 0 && s_e != 0) {
            if (s_e > 0) {
                side = SIDE_AH_ABOVE;
            } else {
                side = SIDE_AH_BELOW;
            }
        } else if  (s0 * s_e == 1) {
            if (s0 > 0) {
                side = SIDE_AH_ABOVE;
            } else {
                side = SIDE_AH_BELOW;
            }
        } else if (s0 * s_e == -1) {
            side = SIDE_AH_RED;
        } else if (s0 != 0 and s_e == 0) {
            if (s0 > 0) {
                side = SIDE_AH_ABOVE;
            } else {
                side = SIDE_AH_BELOW;
            }
        } else {
            MODUS_ASSERT(false);
        }
    } else {
        MODUS_ASSERT(false);
    }
    return side;
}

int sign_vector_to_base3(std::string sv) {
    int n = 0;
    int b = 1;
    for (int i = 0; i < sv.size(); i++) {
        if (sv[i] == '0') {
            n += 0 * b;
        } else if (sv[i] == '+') {
            n += 1 * b;
        } else if (sv[i] == '-') {
            n += 2 * b;
        } else {
            MODUS_ASSERT(false);
        }
        b *= 3;
    }
    return n;
}

int sign_vector_to_base(std::string sv, const Eigen::VectorXi& base) {
    int n = 0;
    int b = 1;
    for (int i = 0; i < sv.size(); i++) {
        if (sv[i] == '-') {
            n += 0 * b;
        } else if (sv[i] == '0') {
            n += 1 * b;
        } else if (sv[i] == '+') {
            n += 2 * b;
        } else {
            MODUS_ASSERT(false);
        }
        b *= base[i];
    }
    return n;
}

bool next_d_face(std::string& sv, const Eigen::VectorXi& sides) {
    int i;
    for (i = 0; i < sv.size(); i++) {
        if (sides[i] == 2) {
            continue;
        }
        if (sv[i] == '-') {
            sv[i] = '+';
            return true;
        } else if (sv[i] == '+') {
            sv[i] = '-';
        }
    }
    return i == sv.size();
}

std::vector<int> initial_hyperplanes(Eigen::MatrixXd& A, 
                                     Eigen::VectorXd& b, 
                                     double eps) {
    // def reorder_halfspaces(A, b, eps=np.finfo(np.float32).eps):
    // A = A.copy()
    // b = b.copy()
    // n = A.shape[0]
    // d = A.shape[1]
    // I = list(range(n))
    // i = 1
    // j = 1
    // while i < n and j < d:
    //     A0 = np.concatenate((A[0:j], A[i,None]), axis=0)
    //     if np.linalg.matrix_rank(A0, eps) > j:
    //         A[j], A[i] = A[i].copy(), A[j].copy()
    //         b[j], b[i] = b[i], b[j]
    //         I[j], I[i] = I[i], I[j]
    //         j += 1
    //     i += 1
    // return A, b, np.array(I, int)
    int n = A.rows();
    int d = A.cols();
    std::vector<int> rindx(n);
    for (int i = 0; i < n; i++) {
        rindx[i] = i;
    }
    int i = 1;
    int j = 1;
    while (i < n && j < d) {
        Eigen::MatrixXd A0(j+1, d);
        A0.topRows(j) = A.topRows(j);
        A0.bottomRows(1) = A.row(i);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
        qr.setThreshold(eps);
        qr.compute(A0);
        if (DEBUG) {
            // std::cout << "i " << i << std::endl;
            // std::cout << "j " << i << std::endl;
            // std::cout << "A0 size " << A0.rows() << " " << A0.cols() << std::endl;
            // std::cout << "A0 rank " << qr.rank() << std::endl;
        }
        if (qr.rank() > j) {
            A.row(j).swap(A.row(i));
            b.row(j).swap(b.row(i));
            int tmp = rindx[j];
            rindx[j] = rindx[i];
            rindx[i] = tmp;
            j += 1;
        }
        i += 1;
    }
    if (DEBUG) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
        qr.setThreshold(eps);
        qr.compute(A.topRows(d));
        MODUS_ASSERT(qr.rank() == d);
    }
    return rindx;
}

IncidenceGraph* initial_arrangement(
    const Eigen::MatrixXd& A, 
    const Eigen::VectorXd& b, 
    double eps) {

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    // MODUS_ASSERT we are given d linearly independent hyperplanes.
    int n = A.rows();
    int d = A.cols();
    MODUS_ASSERT(n == d);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    if (DEBUG) {
        // std::cout << "rank(A) " << qr.rank() << std::endl;
    }
    MODUS_ASSERT(qr.rank() == d);

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "    init: " << n << "x" << n << std::endl;
        std::cout << "    rank: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
    }

    // Allocate k-faces, 0 <= k <= d.
    IncidenceGraph* I = new IncidenceGraph(d);
    I->_nodes.resize((int) pow(3, d));
    for (int i = 0; i < I->_nodes.size(); i++) {
        Node* node = I->_node_pool.construct(-10);
        node->_id = i;
        node->_graph = I;
        I->_nodes[i] = node;
    }
    // Copy halfspaces.
    I->A = A;
    I->b = b;

    // I->_arc_pool.set_next_size(24000000);

    // We build faces from top to bottom, starting with the d+1 face.
    Node* one = I->make_node(d+1);
    one->sign_vector = "{1}";
    I->add_node_to_rank(one);
    // Next we build the regions, i.e. d faces.
    int num_d_faces = 1 << d;
    std::string sv(d, '+');
    for (int i = 0; i < num_d_faces; i++) {
        // Create sign vector string.
        for (int j = 0; j < d; j++) {
            if (i & 1 << j) {
                sv[j] = '+';
            } else {
                sv[j] = '-';
            }
        }
        if (DEBUG) {
            // std::cout << sv << std::endl;
        }
        // Create region.
        int i_f = sign_vector_to_base3(sv);
        Node* f = I->_nodes[i_f];
        f->rank = d;
        f->_key = sv;
        I->add_node_to_rank(f);
        I->add_arc(f, one);
        I->_nodes[i_f] = std::move(f);
    }
    // k faces, 0 <= k <= d-1
    for (int k = d-1; k >= 0; k--) {
        for (Node* g : I->rank(k + 1)) {
            sv = g->_key;
            MODUS_ASSERT(sv.size() == d);
            for (int i = 0; i < d; i++) {
                if (sv[i] == '0') {
                    continue;
                }
                char tmp = sv[i];
                sv[i] = '0';
                Node* f = I->_nodes[sign_vector_to_base3(sv)];
                f->rank = k;
                f->_key = sv;
                if (f->_black_bit == 0) {
                    I->add_node_to_rank(f);
                    f->_black_bit = 1;
                }
                I->add_arc(f, g);
                sv[i] = tmp;
            }
        }
    }
    // -1 face
    Node* zero = I->make_node(-1);
    Node* v = I->rank(0)[0];
    zero->sign_vector = "{0}";
    I->add_arc(zero, v);
    I->add_node_to_rank(zero);

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   faces: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
    }
    
    // Compute interior point for 0 face.
    v->update_interior_point(eps);

    // Compute interior point for 1 faces.
    {
        for (Node* e : I->rank(1)) {
            e->update_interior_point(eps);
        }
    }

    // Compute interior point for k faces.
    for (int k = 2; k < d + 1; k++) {
        for (Node* f : I->rank(k)) {
            f->update_interior_point(eps);
        }
    }

    // Reset black bit.
    for (int i = 0; i < I->_nodes.size(); i++) {
        I->_nodes[i]->_black_bit = 0;
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  int pt: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        std::cout << "   total: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
        << " ms" << std::endl << std::endl;
    }

    return I;
}

void partial_preprocess_hyperplanes(Eigen::MatrixXd& A, 
                                    Eigen::VectorXd& b,
                                    Eigen::VectorXi& sides,
                                    double eps) {
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    Eigen::VectorXi sides_;
    MODUS_ASSERT(A.size() > 0);

    // if (DEBUG) {
    //     Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //     DLOG(INFO) << "\noriginal" << "\n"
    //                << std::fixed << std::setprecision(6) << std::setfill(' ')
    //                << A.format(CleanFmt) << "\n"
    //                << sides.transpose().format(CleanFmt);
    // }

    // Remove duplicate hyperplanes, favoring sides=2 over sides=3.
    std::vector<int> lex_order = lexographic_argsort(A);
    int n_unique = 1;
    for (int i = 1; i < lex_order.size(); i++) {
        int i0 = lex_order[i-1];
        int i1 = lex_order[i];
        if ((A.row(i0) - A.row(i1)).norm() > eps) {
            n_unique += 1;
        }
    }
    A_.resize(n_unique, A.cols());
    b_.resize(n_unique);
    sides_.resize(n_unique);
    A_.row(0) = A.row(lex_order[0]);
    b_[0] = b[lex_order[0]];
    sides_[0] = sides[lex_order[0]];
    int k = 0;
    for (int i = 1; i < lex_order.size(); i++) {
        int i0 = lex_order[i-1];
        int i1 = lex_order[i];
        // BUG Add in check for b.
        if ((A.row(i0) - A.row(i1)).norm() < eps) {
            sides_[k] = std::min(sides_[k], sides[i1]);
        } else {
            A_.row(++k) = A.row(i1);
            b_[k] = b[i1];
            sides_[k] = sides[i1];
        }
    }
    // if (DEBUG) {
    //     Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //     DLOG(INFO) << "\nunique" << "\n"
    //                << std::fixed << std::setprecision(6) << std::setfill(' ')
    //                << A_.format(CleanFmt) << "\n"
    //                << sides_.transpose().format(CleanFmt);
    // }

    // Reorder hyperplanes so that sides=2 is first.
    k = 0;
    for (int i = 0; i < n_unique; i++) {
        if (sides_[i] == 2) {
            A_.row(i).swap(A_.row(k));
            b_.row(i).swap(b_.row(k));
            sides_.row(i).swap(sides_.row(k));
            k++;
        }
    }
    // if (DEBUG) {
    //     Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //     DLOG(INFO) << "\nside order" << "\n"
    //                << std::fixed << std::setprecision(6) << std::setfill(' ')
    //                << A_.format(CleanFmt) << "\n"
    //                << sides_.transpose().format(CleanFmt);
    // }

    // Reorder hyperplanes to be linearly independent in top d rows.
    std::vector<int> li_order = initial_hyperplanes(A_, b_, eps);
    Eigen::VectorXi sides_2(sides_.size());
    for (int i = 0; i < li_order.size(); i++) {
        sides_2[i] = sides_[li_order[i]];
    }
    // if (DEBUG) {
    //     Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //     DLOG(INFO) << "\nlinearly indepdendent" << "\n"
    //                << std::fixed << std::setprecision(6) << std::setfill(' ')
    //                << A_.format(CleanFmt) << "\n"
    //                << sides_2.transpose().format(CleanFmt);
    // }

    // Remove 0 hyperplanes and output.
    int n_nonzero = 0;
    for (int i = 0; i < A_.rows(); i++) {
        if (A_.row(i).norm() < eps) {
            continue;
        }
        n_nonzero++;
    }
    A.resize(n_nonzero, Eigen::NoChange);
    b.resize(n_nonzero);
    sides.resize(n_nonzero);
    k = 0;
    for (int i = 0; i < A_.rows(); i++) {
        if (A_.row(i).norm() < eps) {
            continue;
        }
        A.row(k) = A_.row(i);
        b.row(k) = b_.row(i);
        sides.row(k) = sides_2.row(i);
        k++;
    }
    // if (DEBUG) {
    //     Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //     DLOG(INFO) << "\nnonzero" << "\n"
    //                << std::fixed << std::setprecision(6) << std::setfill(' ')
    //                << A.format(CleanFmt) << "\n"
    //                << sides.transpose().format(CleanFmt);
    // }

    // // Output.
    // A = A_;
    // b = b_;
    // sides = sides_2;
}

IncidenceGraph* partial_initial_arrangement(const Eigen::MatrixXd& A,
                                            const Eigen::VectorXd& b,
                                            const Eigen::VectorXi& sides,
                                            double eps) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    // Assert we are given d linearly independent hyperplanes.
    int n = A.rows();
    int d = A.cols();
    if (DEBUG) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
        qr.setThreshold(eps);
        qr.compute(A.topRows(d));
        // if (qr.rank() != d) {
        //     DLOG(ERROR) << "rank " << qr.rank() << " != " << d;
        // }
        MODUS_ASSERT(qr.rank() == d);
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "    init: " << d << "x" << d << std::endl;
        std::cout << "    rank: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
    }

    // Calculate the number of faces.
    int num_faces = 1;
    for (int i = 0; i < d; i++) {
        num_faces *= sides[i];
        if (DEBUG) {
            MODUS_ASSERT(sides[i] == 2 || sides[i] == 3);
        }
    }

    // Allocate k-faces, 0 <= k <= d.
    IncidenceGraph* I = new IncidenceGraph(d);
    I->_nodes.resize((int) num_faces);
    for (int i = 0; i < I->_nodes.size(); i++) {
        Node* node = I->_node_pool.construct(-10);
        node->_id = i;
        node->_graph = I;
        I->_nodes[i] = node;
    }

    // Copy halfspaces.
    I->A = A.topRows(d);
    I->b = b.topRows(d);

    // We build faces from top to bottom, starting with the d+1 face.
    Node* one = I->make_node(d+1);
    one->sign_vector = "{1}";
    I->add_node_to_rank(one);

    // Calculate the number of d faces.
    int num_d_faces = 1;
    for (int i = 0; i < d; i++) {
        num_d_faces *= sides[i] - 1;
    }

    // Next we build the regions, i.e. d faces.
    std::string sv(d, '-');
    int count = 0;
    for (int i = 0; i < num_d_faces; i++) {
        // Create region.
        int i_f = sign_vector_to_base(sv, sides);
        Node* f = I->_nodes[i_f];
        f->rank = d;
        f->_key = sv;
        I->add_node_to_rank(f);
        I->add_arc(f, one);
        I->_nodes[i_f] = std::move(f);

        // Get next d face.
        bool succ = next_d_face(sv, sides);
        MODUS_ASSERT(succ);

        count++;
    }
    MODUS_ASSERT(count == num_d_faces);

    // Build k faces, 0 <= k <= d-1
    for (int k = d-1; k >= 0; k--) {
        for (Node* g : I->rank(k + 1)) {
            sv = g->_key;
            MODUS_ASSERT(sv.size() == d);
            for (int i = 0; i < d; i++) {
                if (sv[i] == '0') {
                    continue;
                }
                char tmp = sv[i];
                sv[i] = '0';
                Node* f = I->_nodes[sign_vector_to_base(sv, sides)];
                f->rank = k;
                f->_key = sv;
                if (f->_black_bit == 0) {
                    I->add_node_to_rank(f);
                    f->_black_bit = 1;
                    count ++;
                }
                I->add_arc(f, g);
                sv[i] = tmp;
            }
        }
    }
    MODUS_ASSERT(count == num_faces);

    // -1 face
    Node* zero = I->make_node(-1);
    zero->sign_vector = "{0}";
    Node* v = I->rank(0)[0];
    I->add_arc(zero, v);
    I->add_node_to_rank(zero);

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   faces: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
    }

    // if (DEBUG) {
    //     for (Node* u : I->_nodes) {
    //         std::cout << u->_key << std::endl;
    //     }
    // }

    // Compute interior point for 0 face.
    v->update_interior_point(eps);

    // Compute interior point for 1 faces.
    {
        for (Node* e : I->rank(1)) {
            e->update_interior_point(eps);
        }
    }

    // Compute interior point for k faces.
    for (int k = 2; k < d + 1; k++) {
        for (Node* f : I->rank(k)) {
            f->update_interior_point(eps);
        }
    }
    I->update_sign_vectors(eps);

    // Reset black bit.
    for (int i = 0; i < I->_nodes.size(); i++) {
        I->_nodes[i]->_black_bit = 0;
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  int pt: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;

        std::cout << "   total: " << 
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
        << " ms" << std::endl << std::endl;
    }

    return I;
}

IncidenceGraph* partial_initial_convex(const Eigen::MatrixXd& A,
                                       const Eigen::VectorXd& b,
                                       double eps) {
    // 1. Find a strictly interior point (if it exists).
    // 2. Partition hyperplanes.
    // 3. Project into null space, affine space, and dual space.
    // 4. Build convex hull and face lattice.
    // 5. Update interior points.

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        std::cout << "    init: " << A.rows() << "x" << A.cols() << std::endl;
        start = std::chrono::high_resolution_clock::now();
        start_total = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 1. Find a strictly interior point.                                     *
     **************************************************************************/
    // Get dimensions.
    int n = A.rows();
    int d = A.cols();
    // Find a strictly interior point.
    Eigen::VectorXd r = interior_point(A, b, eps);
    if (DEBUG) {
        // DLOG(INFO) << "\n\t" << "interior point: " << r.transpose();
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
    Eigen::MatrixXd A_c, A_s, A_proj;
    Eigen::VectorXd b_c, b_s;
    Eigen::VectorXi idx_c, idx_s;
    // Find hyperplanes which are always in contact.
    Eigen::VectorXd res = A*r - b;
    arg_where_zero(res, eps, idx_c);
    arg_where_nonzero(res, eps, idx_s);
    // Partition hyperplanes.
    GetRows(A, idx_c, A_c);
    GetRows(A, idx_s, A_s);
    GetRows(b, idx_c, b_c);
    GetRows(b, idx_s, b_s);

    /**************************************************************************
     * 3. Project into null space, affine space, and dual affine space.       *
     **************************************************************************/
    // Get nullspace. 
    Eigen::MatrixXd kernel;
    kernel_basis(A_c, kernel, eps);
    std::cout << "kernel\n" << kernel << std::endl;
    // Project hyperplanes and interior point into contacting nullspace.
    A_s = A_s * kernel;
    A_c = A_c * kernel;
    A_proj = A * kernel;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.setThreshold(eps);
    qr.compute(kernel);
    Eigen::VectorXd tmp = qr.solve(r); // Project interior point.
    r = tmp;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "\t\t[", "]");
    // DLOG(INFO) << "\n" << "A_proj" << std::fixed << std::setprecision(6) << std::setfill(' ')
    //  << "\nA_proj\n" << A_proj.format(CleanFmt) << std::endl;

    // Get affine subspace.
    Eigen::MatrixXd affine = orth(A_s.transpose(), eps);
    // Project hyperplanes and interior point into affine subspace.
    A_s = A_s * affine;
    A_c = A_c * affine;
    A_proj = A_proj * affine;
    r = affine.transpose() * r;
    
    // DLOG(INFO) << "\n" << "A_proj" << std::fixed << std::setprecision(6) << std::setfill(' ')
    //  << "\nA_proj\n" << A_proj.format(CleanFmt) << "\n" << b.transpose() << std::endl;

    // Get dual points.
    Eigen::VectorXd b_off = b_s - A_s * r;
    Eigen::MatrixXd dual = A_s;
    for (int i = 0; i < A_s.rows(); i++) {
        dual.row(i) /= b_off[i];
    }
    
    // Project dual points into affine subspace.
    Eigen::MatrixXd dual_aff;
    project_affine(dual.transpose(), eps, dual_aff);
    int d_eff = dual_aff.rows();
    std::cout << "dual_aff\n" << std::fixed << std::setprecision(6) << std::setfill(' ')
    << dual_aff.format(CleanFmt) << std::endl;

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "partproj: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }

    /**************************************************************************
     * 4. Build convex hull and face lattice.                                 *
     **************************************************************************/
    // Convex hull.
    Eigen::MatrixXi M = convex_hull(dual_aff, eps);
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "convhull: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    // Build face lattice.
    IncidenceGraph* graph = build_graph_cone(M.transpose(), n, d_eff, eps);
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
    graph->set_hyperplanes(A_s, b_s);
    // graph->set_hyperplanes(A_proj, b);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->update_interior_point(eps);
        }
    }

    // Lift interior points and update sign vectors.
    graph->set_hyperplanes(A, b);
    for (int k = 0; k <= d_eff + 1; k++) {
        for (Node* u : graph->rank(k)) {
            u->interior_point = kernel * affine * u->interior_point;
        }
    }
    graph->update_sign_vectors(eps);

    std::vector<std::string> sv = graph->get_sign_vectors();
    std::cout << "convex" << std::endl;
    for (std::string s : sv) {
        std::cout << s << std::endl;
    }

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "  update: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_total).count() / 1e6
        << " ms" << std::endl << std::endl;
    }

    return graph;
}

void increment_arrangement(Eigen::VectorXd a, double b, 
                           IncidenceGraph* I, double eps) {
    // Normalize halfspace, |a| = 1.
    double norm_a = a.norm();
    a /= norm_a;
    b /= norm_a;
    I->add_hyperplane(a, b);
    I->update_sign_vectors(eps);

    // =========================================================================
    // Phase 1: Find an edge e₀ in 𝓐(H) such that cl(e₀) ∩ h ≠ ∅
    // =========================================================================

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_total;
    if (PROFILE) {
        std::cout << "    incr: " << I->A.rows() << "x" << I->A.cols() << std::endl;
        std::cout << " phase 1: " << std::endl;
        start = std::chrono::high_resolution_clock::now();
    }
    int n = I->A.rows();
    Node* v = I->rank(0)[0];
    // Find an incident edge e on v such that aff(e) is not parallel to h.
    Node* e;
    Eigen::VectorXd v_e;
    double dist;
    for (auto it = v->superfaces.begin(); it != v->superfaces.end(); it++) {
        e = *it;
        get_vector(e, v_e);
        dist = std::abs(a.dot(v_e));
        if (dist > eps) {
            break;
        }
    }
    if (DEBUG) {
        // std::cout << "\tdot: " << dist << std::endl;
        MODUS_ASSERT(dist > eps);
    }
    // Find edge e₀ such that cl(e₀) ∩ h ≠ ∅.
    Node* e0 = e;
    Eigen::VectorXd v_e0 = v_e;
    v_e0.normalize();
    while (true) {
        if (get_color_edge(e0, a, b, eps) > COLOR_AH_WHITE) {
            break;
        }
        // Find v(e0) closer to h.
        Node* v;
        if (e0->subfaces.size() == 2) {
            auto i = e0->subfaces.begin();
            Node* v0 = *i++;
            Node* v1 = *i;
            double d0 = abs(a.dot(v0->interior_point) - b);
            double d1 = abs(a.dot(v1->interior_point) - b);
            if (d0 < d1) {
                v = v0;
            } else {
                v = v1;
            }
        } else if (e0->subfaces.size() == 1) {
            v = *e0->subfaces.begin();
        } else {
            MODUS_ASSERT(false);
        }
        // Find e in v such that aff(e0) == aff(e).
        Node* e_min;
        Eigen::VectorXd v_min;
        double min_dist = std::numeric_limits<float>::infinity();
        if (DEBUG) {
            // std::cout << "v.superfaces: " << v->superfaces.size() << std::endl;
        }
        for (Node* e : v->superfaces) {
            if (e == e0) {
                continue;
            }
            get_vector(e, v_e);
            v_e.normalize();
            double dist = (v_e0 + v_e).norm();
            // double dist = (v_e - v_e0.dot(v_e) * v_e).norm();
            if (dist < min_dist) {
                e_min = e;
                v_min = v_e;
                min_dist = dist;
            }
            if (DEBUG) {
                // std::cout << e->_key << std::endl;
                // std::cout << (v_e0 + v_e).norm() << std::endl;
                // std::cout << dist << std::endl;
            }
        }
        e0 = e_min;
        if (DEBUG) {
            // std::cout << " e min: " << e_min->_key << std::endl;
            // std::cout << " d min: " << min_dist << std::endl;
        }
    }
    if (DEBUG) {
        // std::cout << "\te0.color: " << 
        //     get_color_ah_string(get_color_edge(e0, a, b, eps)) << std::endl;
        // std::cout << "\te0.key  : " << e0->_key << std::endl;
    }

    // =========================================================================
    // Phase 2: Mark all faces f with cl(f) ∩ h ≠ ∅ pink, red, or crimson.
    // =========================================================================

    int mark_count = 0;
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::cout << " phase 2: " << std::endl;
    }
    // Add some 2 face incident upon e₀ to Q and mark it green.
    Node* f = *e0->superfaces.begin();
    f->_color = COLOR_AH_GREEN;
    std::list<Node*> Q;
    Q.push_back(f);
    // Color vertices, edges, and 2 faces of 𝓐(H).
    int d = a.size();
    std::vector<std::vector<Node*> > L(d+1);
    for (int k = 0; k < d+1; k++) {
        L[k].clear();
    }
    // Node* g;
    while (!Q.empty()) {
        f = Q.front();
        Q.pop_front();
        for (Node* e : f->subfaces) {
            // Skip edges we have already seen.
            if (e->_color != COLOR_AH_WHITE) {
                continue;
            }
            int color_e = get_color_edge(e, a, b, eps);
            if (color_e > COLOR_AH_WHITE) {
                // Mark each white vertex v ∈ h crimson and insert v into L₀.
                for (Node* v : e->subfaces) {
                    if (v->_color == COLOR_AH_WHITE) {
                        int color_v = get_color_vertex(v, a, b, eps);
                        if (color_v == COLOR_AH_CRIMSON) {
                            v->_color = color_v;
                            L[0].push_back(v);
                        }
                    }
                }
                // Color e and insert e into L₁.
                e->_color = color_e;
                L[1].push_back(e);
                // Mark all white 2 faces green and put them into Q.
                for (Node* g : e->superfaces) {
                    if (g->_color == COLOR_AH_WHITE) {
                        g->_color = COLOR_AH_GREEN;
                        Q.push_back(g);
                    }
                }
                // Mark the side of the hyperplane the edge is on.
                e->_side = get_side_edge(e, a, b, eps);
            }
        }
    }
    // Color k faces, 2 ≤ k ≤ d.
    for (int k = 2; k < d+1; k++) {
        for (Node* f : L[k-1]) {
            for (Node* g : f->superfaces) {
                if (g->_color != COLOR_AH_WHITE && g->_color != COLOR_AH_GREEN) {
                    continue;
                }
                // Case 1: f is pink. Mark g red if it has a red subface or pink
                // subfaces on both sides of h. Otherwise, mark g pink.
                if (f->_color == COLOR_AH_PINK) {
                    bool above = false;
                    bool below = false;
                    for (Node* f_g : g->subfaces) {
                        if (f_g->_color == COLOR_AH_CRIMSON || f_g->_color == COLOR_AH_WHITE) {
                            continue;
                        }
                        if (f_g->_color == COLOR_AH_RED) {
                            above = true;
                            below = true;
                            break;
                        }
                        MODUS_ASSERT(f_g->_color == COLOR_AH_PINK);
                        MODUS_ASSERT(f_g->_side == SIDE_AH_ABOVE || f_g->_side == SIDE_AH_BELOW);
                        if (f_g->_side == SIDE_AH_ABOVE) {
                            above = true;
                        } else if (f_g->_side == SIDE_AH_BELOW) {
                            below = true;
                        }
                        if (DEBUG) {
                            int sign = get_position(a.dot(f_g->interior_point) - b, eps);
                            int side = f_g->_side;
                            MODUS_ASSERT_WARNING(sign != 0,
                                "Warning: pink face %d is possibly red, check numerical stability",
                                f_g->_id);
                            MODUS_ASSERT_WARNING(sign > 0 ? side == SIDE_AH_ABOVE : side == SIDE_AH_BELOW,
                                "Warning: pink face %d is possibly red, check numerical stability",
                                f_g->_id);
                        }
                    }
                    if (above && below) {
                        g->_color = COLOR_AH_RED;
                    } else {
                        g->_color = COLOR_AH_PINK;
                        g->_side = above ? SIDE_AH_ABOVE : SIDE_AH_BELOW;
                        MODUS_ASSERT(above != below);
                    }
                }
                else if (f->_color == COLOR_AH_RED) {
                    g->_color = COLOR_AH_RED;
                }
                else if (f->_color == COLOR_AH_CRIMSON) {
                    bool above = false;
                    bool below = false;
                    bool crimson = true;
                    MODUS_ASSERT(g->subfaces.size() > 0);
                    for (Node* f_g : g->subfaces) {
                        if (f_g->_color != COLOR_AH_CRIMSON) {
                            crimson = false;
                        }
                        // For pink subfaces, use precomputed side.
                        if (f_g->_color == COLOR_AH_PINK) {
                            above |= f_g->_side == SIDE_AH_ABOVE;
                            below |= f_g->_side == SIDE_AH_BELOW;
                            MODUS_ASSERT(f_g->_side == SIDE_AH_ABOVE || f_g->_side == SIDE_AH_BELOW);
                        }
                        // For white subfaces, compute the side.
                        if (f_g->_color == COLOR_AH_WHITE) {
                            int s = get_position(a.dot(f_g->interior_point) - b, eps);
                            MODUS_ASSERT(s != 0);
                            above |= (s ==  1);
                            below |= (s == -1);
                        }
                    }
                    if (crimson) {
                        g->_color = COLOR_AH_CRIMSON;
                    } else {
                        g->_color = COLOR_AH_PINK;
                        g->_side = above ? SIDE_AH_ABOVE : SIDE_AH_BELOW;
                        if (above == below) {
                            std::cout << "A\n" << I->A << std::endl;
                        }
                        MODUS_ASSERT(above != below,
                            "Error: pink face %d with crimson subface and subfaces %s%s%s, check numerical stability",
                            g->_id, (above ? "above" : ""), (above && below ? " and " : ""), (below ? "below" : ""));
                    }
                }
                else {
                    std::cout << get_color_ah_string(f->_color) << std::endl;
                    MODUS_ASSERT(false);
                }
                L[k].push_back(g);
            }
        }
    }

    // Verify face coloring satisfies Observation 7.2 and Table 7.3 of
    // Edelsbrunner.
    if (DEBUG)
    {
        // 7.2 (i) Vertices are either white or crimson.
        for (Node* v : I->rank(0)) {
            MODUS_ASSERT(v->_color == COLOR_AH_WHITE || v->_color == COLOR_AH_CRIMSON,
                "Error: vertex %d has color %s", v->_id, get_color_ah_string(v->_color).c_str());
        }
        for (size_t k = 1; k < d + 1; k++) {
            for (Node* g : I->rank(k)) {
                if (g->_color == COLOR_AH_WHITE || g->_color == COLOR_AH_GREEN) {
                    // Table 7.3 [:,1] The subfaces of a white face are all
                    // white.
                    for (Node* f : g->subfaces) {
                        MODUS_ASSERT(f->_color == COLOR_AH_WHITE,
                            "Error: white superface %d has subface %d with color %s", 
                            g->_id, f->_id, get_color_ah_string(f->_color).c_str());
                    }
                    // Table 7.3 [1,:] The superfaces of a white face are white,
                    // pink, or red.
                    if (k == d) continue;
                    for (Node* h : g->superfaces) {
                        MODUS_ASSERT(h->_color == COLOR_AH_WHITE || h->_color == COLOR_AH_PINK || h->_color == COLOR_AH_RED,
                            "Error: white subface %d has crimson superface %d", g->_id, h->_id);
                    }
                }
                else if (g->_color == COLOR_AH_PINK) {
                    // Table 7.3 [:,2] The subfaces of a pink face are either
                    // white, pink, or crimson.
                    for (Node* f : g->subfaces) {
                        MODUS_ASSERT(f->_color == COLOR_AH_WHITE || f->_color == COLOR_AH_PINK || f->_color == COLOR_AH_CRIMSON,
                            "Error: pink superface %d has subface %d with color %s",
                            g->rank, f->rank, get_color_ah_string(f->_color).c_str());
                    }
                    // Table 7.3 [2,:] The superfaces of a pink face are either
                    // pink or red.
                    if (k == d) continue;
                    for (Node* h : g->superfaces) {
                        MODUS_ASSERT(h->_color == COLOR_AH_PINK || h->_color == COLOR_AH_RED,
                            "Error: pink subface %d has superface %d with color %s",
                            g->rank, h->rank, get_color_ah_string(h->_color).c_str());
                    }
                }
                else if (g->_color == COLOR_AH_RED) {
                    // Table 7.3 [:,3] The subfaces of a red face are either
                    // white, pink, or red.
                    for (Node* f : g->subfaces) {
                        MODUS_ASSERT(f->_color == COLOR_AH_WHITE || f->_color == COLOR_AH_PINK || f->_color == COLOR_AH_RED,
                            "Error: red superface %d has subface %d with color %s",
                            g->rank, f->rank, get_color_ah_string(f->_color).c_str());
                    }
                    // Table 7.3 [3,:] The superfaces of a red face are red.
                    if (k == d) continue;
                    for (Node* h : g->superfaces) {
                        MODUS_ASSERT(h->_color == COLOR_AH_RED,
                            "Error: red subface (%d) rank %d/%d has superface (%d) rank %d with color %s",
                            g->_id, g->rank, d, h->_id, h->rank, get_color_ah_string(h->_color).c_str());
                    }
                }
                else if (g->_color == COLOR_AH_CRIMSON) {
                    // Table 7.3 [:,4] The subfaces of a crimson face are
                    // crimson.
                    for (Node* f : g->subfaces) {
                        MODUS_ASSERT(f->_color == COLOR_AH_CRIMSON,
                            "Error: crimson superface %d has subface %d with color %s",
                            g->rank, f->rank, get_color_ah_string(f->_color).c_str());
                    }
                    // Table 7.3 [4,:] The superfaces of a crimson face are pink
                    // or crimson.
                    if (k == d) continue;
                    for (Node* h : g->superfaces) {
                        MODUS_ASSERT(h->_color == COLOR_AH_PINK || h->_color == COLOR_AH_CRIMSON,
                            "Error: crimson subface (%d) rank %d has superface (%d) rank %d with color %s",
                            g->_id, g->rank, h->_id, h->rank, get_color_ah_string(h->_color).c_str());
                    }
                }
            }
        }
    }

    // =========================================================================
    // Phase 3: Update all marked faces.
    // =========================================================================

    double time_step_1 = 0;
    double time_step_2 = 0;
    double time_step_3 = 0, time_step_3_copy = 0, time_step_3_ins = 0;
    double time_step_4 = 0;
    double time_step_5 = 0;
    double time_step_6 = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_step;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_1;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_2;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_3;
    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::cout << " phase 3:" << std::endl;
    }
    // I->_arc_pool.set_next_size(24000000);

    for (int k = 0; k < d + 1; k++) {
        int n_Lk = L[k].size();
        for (int i = 0; i < n_Lk; i++) {
            Node* g = L[k][i];
            switch (g->_color) {

                case COLOR_AH_PINK: {
                g->_color = COLOR_AH_GREY;
                // Add to grey subfaces of superfaces.
                for (Node* u : g->superfaces) {
                    u->_grey_subfaces.push_back(g);
                }
                break; }

                case COLOR_AH_CRIMSON: {
                g->_color = COLOR_AH_BLACK;
                // Add to black subfaces of superfaces.
                for (Node* u : g->superfaces) {
                    u->_black_subfaces.push_back(g);
                }
                break; }
                
                case COLOR_AH_RED: {

                if (PROFILE) {
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 1. Create g_a = g ∩ h⁺ and g_b = g ∩ h⁻. Remove g from
                // 𝓐(H) and Lₖ and replace with g_a, g_b.
                if (DEBUG) {
                    // std::cout << "Splitting the following node..." << std::endl;
                    // std::cout << *g << std::endl;
                }
                Node* g_a = g;
                g_a->_color = COLOR_AH_GREY;
                Node* g_b = I->make_node(k);
                g_b->_color = COLOR_AH_GREY;
                I->add_node_to_rank(g_b);
                L[k].push_back(g_b);

                g->update_sign_vector(eps);
                g_a->_key = g->sign_vector;
                g_a->_key.back() = '+';
                g_a->_side = SIDE_AH_ABOVE;
                g_b->_key = g->sign_vector;
                g_b->_key.back() = '-';
                g_b->_side = SIDE_AH_BELOW;

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_1 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 2. Create the black face f = g ∩ h, connect it to g_a
                // and g_b, and put f into 𝓐(H) and Lₖ₋₁.
                Node* f = I->make_node(k-1);
                f->_color = COLOR_AH_BLACK;
                I->add_arc(f, g_a);
                I->add_arc(f, g_b);
                g_a->_black_subfaces = {f};
                g_b->_black_subfaces = {f};
                L[k-1].push_back(f);
                I->add_node_to_rank(f);

                f->_key = g->sign_vector;
                f->_key.back() = '0';

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_2 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 3. Connect each red superface of g with g_a and g_b.
                for (Node* r : g->superfaces) {
                    if (DEBUG) {
                        MODUS_ASSERT(r->_color == COLOR_AH_RED || r->rank == d+1);
                    }
                    I->add_arc(g_b, r);
                    r->_grey_subfaces.push_back(g_a);
                    r->_grey_subfaces.push_back(g_b);
                }

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_3 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 4. Connect each white or grey subface of g with g_a if
                // it is in h⁺, and with g_b, otherwise.
                ArcListIterator iter = g->subfaces.begin();
                ArcListIterator end  = g->subfaces.end();
                while (iter != end) {
                    Arc* arc = iter.arc;
                    Node* u = *iter++;
                    if (DEBUG) {
                        // std::cout << i_u << std::endl;
                        // std::cout << *iter.arc << std::endl;
                    }

                    if (u->_color == COLOR_AH_BLACK) {
                        continue;
                    }
                    // if (u->_color != COLOR_AH_WHITE && u->_color != COLOR_AH_GREY) {
                    //     if (DEBUG) {
                    //         MODUS_ASSERT(u->_color == COLOR_AH_BLACK);
                    //     }
                    //     continue;
                    // }
                    int side = 0;
                    if (u->_color == COLOR_AH_GREY) {
                        if (u->_side == SIDE_AH_UNSET) {
                            int sign = get_position(a.dot(u->interior_point) - b, eps);
                            if (sign == 1) {
                                side = SIDE_AH_ABOVE;
                            } else if (sign == -1) {
                                side = SIDE_AH_BELOW;
                            } else {
                                MODUS_ASSERT(sign != 0,
                                    "Error: grey face is actually black, check numerical stability");
                            }
                        } else {
                            MODUS_ASSERT(u->_side == SIDE_AH_ABOVE || u->_side == SIDE_AH_BELOW);
                            side = u->_side;
                        }
                    } else if (u->_color == COLOR_AH_WHITE) {
                        int sign = get_position(a.dot(u->interior_point) - b, eps);
                        if (sign == 1) {
                            side = SIDE_AH_ABOVE;
                        } else if (sign == -1) {
                            side = SIDE_AH_BELOW;
                        } else {
                            MODUS_ASSERT(sign != 0,
                                    "Error: white face is actually black, check numerical stability");
                        }
                    }
                    MODUS_ASSERT(side == SIDE_AH_ABOVE || side == SIDE_AH_BELOW);
                    if (side == SIDE_AH_ABOVE) {
                        if (u->_color == COLOR_AH_GREY) {
                            g_a->_grey_subfaces.push_back(u);
                        }
                    } else if (side == SIDE_AH_BELOW) {
                        if (u->_color == COLOR_AH_GREY) {
                            g_b->_grey_subfaces.push_back(u);
                        }
                        g->subfaces._remove_arc(arc);
                        u->superfaces._remove_arc(arc->_dst_arc);
                        I->add_arc(u, g_b, arc, arc->_dst_arc);
                    }
                }

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_4 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 5. If k = 1, connect f with the -1 face, and connect f
                // with the black subfaces of the grey subfaces of g, otherwise.
                if (k == 1) {
                    Node* zero = I->rank(-1)[0];
                    I->add_arc(zero, f);
                } else {
                    // std::set<Node*> V;
                    std::vector<Node*> V;
                    V.clear();
                    for (Node* u : g->_grey_subfaces) {
                        for (Node* v : u->_black_subfaces) {
                            if (PROFILE) {
                                // TODO
                            }
                            if (v->_black_bit == 0) {
                                V.push_back(v);
                                v->_black_bit = 1;
                            }
                        }
                    }
                    if (PROFILE) {
                        // TODO
                    }
                    for (Node* v : V) {
                        v->_black_bit = 0;
                        I->add_arc(v, f);
                    }
                }

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_5 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }
                // Step 6. Update the interior points for f, g_a, and g_b.
                f->update_interior_point(eps);
                if (DEBUG) {
                    
                }
                g_a->update_interior_point(eps);
                if (DEBUG) {
                    
                }
                g_b->update_interior_point(eps);
                if (DEBUG) {
                    // std::cout << "f" << std::endl;
                    // std::cout << *f << std::endl;
                    // std::cout << "g_a" << std::endl;
                    // std::cout << *g_a << std::endl;
                    // std::cout << "g_b" << std::endl;
                    // std::cout << *g_b << std::endl;
                }

                if (PROFILE) {
                    auto end_step = std::chrono::high_resolution_clock::now();
                    time_step_6 += std::chrono::duration_cast<std::chrono::nanoseconds>(end_step - start_step).count() / 1e6;
                    start_step = std::chrono::high_resolution_clock::now();
                }

                break; }
                default: {
                    std::cout << "g: " << get_color_ah_string(g->_color) << std::endl;
                    MODUS_ASSERT(false);
                }
            }
        }
    }

    for (int k = 0; k < d + 1; k++) {
        for (Node* u : L[k]) {
            u->_color = COLOR_AH_WHITE;
            u->_side = SIDE_AH_UNSET;
            u->_grey_subfaces.clear();
            u->_black_subfaces.clear();
        }
    }
    I->update_sign_vectors(eps);

    if (PROFILE) {
        auto end = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < d + 1; k++) {
            std::cout << "\tL(" << k << "):" << L[k].size() << std::endl;
        }
        std::cout << "  step 1: " << time_step_1 << " ms" << std::endl;
        std::cout << "  step 2: " << time_step_2 << " ms" << std::endl;
        std::cout << "  step 3: " << time_step_3 << " ms" << std::endl;
        // std::cout << " copy 3: " << time_step_3_copy << " ms" << std::endl;
        // std::cout << "  ins 3: " << time_step_3_ins << " ms" << std::endl;
        std::cout << "  step 4: " << time_step_4 << " ms" << std::endl;
        std::cout << "  step 5: " << time_step_5 << " ms" << std::endl;
        std::cout << "  step 6: " << time_step_6 << " ms" << std::endl;
        std::cout << "   total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e6
        << " ms" << std::endl << std::endl;
    }
    
}