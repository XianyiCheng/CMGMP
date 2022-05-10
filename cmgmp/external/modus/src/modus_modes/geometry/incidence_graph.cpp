#include <modus/common/eigen.hpp>
#include <modus/modes/geometry/incidence_graph.hpp>
#include <modus/common/linear_algebra.hpp>
#include <modus/common/assert.hpp>
#include <iostream>
#include <iomanip>
// #include <glog/logging.h>


static int DEBUG=0;
static int WARN_TOL=0;

using namespace modus;

void print_vector_int(const std::vector<int>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
}

std::ostream& operator<<(std::ostream& out, const std::vector<int>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        out << vec[i] << " ";
    }
    return out;
}

std::string get_color_ah_string(int color) {
    switch (color) {
        case COLOR_AH_WHITE:    return "WHITE";
        case COLOR_AH_PINK:     return "PINK";
        case COLOR_AH_RED:      return "RED";
        case COLOR_AH_CRIMSON:  return "CRIMSON";
        case COLOR_AH_GREY:     return "GREY";
        case COLOR_AH_BLACK:    return "BLACK";
        case COLOR_AH_GREEN:    return "GREEN";
        default:                return "INVALID";
    }
}

int get_position(double v, double eps) {
    MODUS_ASSERT(eps > 0);
    // if (WARN_TOL && std::abs(v) < 2 * eps && 0.5 * eps < std::abs(v)) {
    //     // LOG(WARNING) << "Approaching tolerance limits: " 
    //     // << 0.5*eps << " < " << std::abs(v) << " < " << 2*eps << std::endl;
    // }
    if (v > eps) {
        return 1;
    } else if (v < -eps) {
        return -1;
    } else {
        return 0;
    }
}

char get_sign(double v, double eps) {
    MODUS_ASSERT(eps > 0);
    // if (WARN_TOL && std::abs(v) < 2 * eps && 0.5 * eps < std::abs(v)) {
    //     LOG(WARNING) << "Approaching tolerance limits: " 
    //     << 0.5*eps << " < " << std::abs(v) << " < " << 2*eps << std::endl;
    // }
    if (v > eps) {
        return '+';
    } else if (v < -eps) {
        return '-';
    } else {
        return '0';
    }
}

void get_position(const Eigen::VectorXd& v, Eigen::VectorXi& pos, double eps) {
    eps = std::abs(eps);
    pos.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
        pos[i] = get_position(v[i], eps);
    }
}

Eigen::VectorXi get_position(const Eigen::VectorXd& v, double eps) {
    Eigen::VectorXi position;
    get_position(v, position, eps);
    return position;
}

void get_sign_vector(const Eigen::VectorXd& v, std::string& sv, double eps) {
    eps = abs(eps);
    sv.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
        sv[i] = get_sign(v[i], eps);
    }
}

std::string get_sign_vector(const Eigen::VectorXd& v, double eps) {
    std::string sv;
    get_sign_vector(v, sv, eps);
    return sv;
}

void arg_where(const std::string& sv, char s, Eigen::VectorXi& idx) {
    int cnt = 0;
    for (int i = 0; i < sv.size(); i++) {
        cnt += sv[i] == s;
    }
    idx.resize(cnt);
    int k = 0;
    for (int i = 0; i < sv.size(); i++) {
        if (sv[i] == s) {
            idx[k++] = i;
        }
    }
}

void arg_equal(const std::string& a, 
               const std::string& b, 
               Eigen::VectorXi& idx) 
{
    MODUS_ASSERT(a.size() == b.size());
    int cnt = 0;
    for (int i = 0; i < a.size(); i++) {
        cnt += a[i] == b[i];
    }
    idx.resize(cnt);
    int k = 0;
    for (int i = 0; i < a.size(); i++) {
        if (a[i] == b[i]) {
            idx[k++] = i;
        }
    }
}

void arg_not_equal(const std::string& a, 
                   const std::string& b, 
                   Eigen::VectorXi& idx)
{
    MODUS_ASSERT(a.size() == b.size());
    int cnt = 0;
    for (int i = 0; i < a.size(); i++) {
        cnt += a[i] != b[i];
    }
    idx.resize(cnt);
    int k = 0;
    for (int i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            idx[k++] = i;
        }
    }
}

int partial_order(const std::string& lhs, const std::string& rhs) {
    if (lhs.size() != rhs.size()) {
        return INCOMPARABLE;
    }
    int less = 0;
    int greater = 0;
    int equal = 0;
    for (int i = 0; i < lhs.size(); i++) {
        if (lhs[i] == '+' && rhs[i] == '-') {
            return INCOMPARABLE;
        }
        else if (lhs[i] == '-' && rhs[i] == '+') {
            return INCOMPARABLE;
        }
        else if (lhs[i] == '0') {
            if (rhs[i] == '0') {
                equal += 1;
            }
            else {
                less += 1;
            }
        }
        else if (rhs[i] == '0') {
            if (lhs[i] == '0') {
                equal += 1;
            } else {
                greater += 1;
            }
        }
        else {
            equal += 1;
        }
        if (less > 0 && greater > 0) {
            return INCOMPARABLE;
        }
    }
    if (less > 0) {
        return STRICTLY_LESS;
    }
    else if (greater > 0) {
        return STRICTLY_GREATER;
    } else {
        return EQUAL;
    }
}

std::vector<int> intersect_sorted(const std::vector<int>& a, 
                                  const std::vector<int>& b) {
    // 
    int n = a.size();
    int m = b.size();
    int i = 0; 
    int j = 0;
    std::vector<int> c; c.clear();
    while (i < n && j < m) {
        if (a[i] < b[j]) {
            i += 1;
        } else if (b[j] < a[i]) {
            j += 1;
        } else {
            c.push_back(a[i]);
            i += 1;
            j += 1;
        }
    }
    return c;
}

std::vector<int> difference_sorted(const std::vector<int>& a,
                                   const std::vector<int>& b) {
    int n = a.size();
    int m = b.size();
    int i = 0;
    int j = 0;
    std::vector<int> c; c.clear();
    while (i < n) {
        if (j == m) {
            c.push_back(a[i]);
            i += 1;
        } else if (b[j] < a[i]) {
            j += 1;
        } else if (a[i] < b[j]) {
            c.push_back(a[i]);
            i += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    return c;
}

std::vector<int> insert_sorted(const std::vector<int>& h, int v) {
    // std::cout << "h ";
    // print_vector_int(h);
    // std::cout << std::endl;
    std::vector<int> a; a.clear();
    bool pushed = false;
    for (int i = 0; i < h.size(); i++) {
        if (v < h[i] && !pushed) {
            a.push_back(v);
            pushed = true;
        }
        a.push_back(h[i]);
    }
    // If the array is empty or v is > all h.
    if (!pushed) {
        a.push_back(v);
    }
    // std::cout << a.size() << std::endl;
    if (DEBUG) {
        MODUS_ASSERT(a.size() == h.size() + 1);
        for (int i = 0; i + 1 < a.size(); i++) {
            if (a[i] > a[i+1]) {
                // std::cout << a[i] << std::endl;
                // std::cout << a[i+1] << std::endl;
            }
            MODUS_ASSERT(a[i] <= a[i+1]);
        }        
    }
    return a;
}

std::vector<int> intersect_preimages(const std::vector<int>& A, 
                                     const std::vector<std::vector<int> >& B) {
    if (A.size() == 0) {
        return std::vector<int>();
    }
    std::vector<int> I = B[A[0]];
    for (int i = 1; i < A.size(); i++) {
        I = intersect_sorted(I, B[A[i]]);
    }
    return I;
}

std::vector<int> closure(const std::vector<int>& S, 
                         const std::vector<std::vector<int> >& V,
                         const std::vector<std::vector<int> >& F) {
    return intersect_preimages(intersect_preimages(S, F), V);
}

IncidenceGraph* build_graph_cone(const Eigen::MatrixXi& M, int n, int d, double eps) {
    // We build the incidence graph from the polar of M. Therefore, columns are
    // vertices and rows are facets.
    int n_verts = M.cols();
    int n_facets = M.rows();

    // Convert M to sorted sparse format.
    std::vector<std::vector<int> > F;   // F(v) -> {f}
    std::vector<std::vector<int> > V;   // V(f) -> {v}
    F.clear(); V.clear();
    for (int i = 0; i < n_verts; i++) {
        std::vector<int> f; f.clear();
        for (int j = 0; j < n_facets; j++) {
            if (M(j,i)) {
                f.push_back(j);
            }
        }
        F.push_back(f);
    }
    for (int i = 0; i < n_facets; i++) {
        std::vector<int> v; v.clear();
        for (int j = 0; j < n_verts; j++) {
            if (M(i,j)) {
                v.push_back(j);
            }
        }
        V.push_back(v);
    }

    // Create incidence graph.
    IncidenceGraph* I = new IncidenceGraph(d+1);

    // Create d+2 rank node, i.e. {1}.
    Node* one = I->make_node(d+2);
    one->_key = std::string();
    one->_vertices.clear();
    one->sign_vector = "{1}";
    I->add_node_to_rank(one);

    // Create d+1 rank node, i.e. the interior of the polytope.
    Node* interior = I->make_node(d+1);
    interior->_key = std::string(n, '-');
    interior->_vertices.clear();
    I->add_node_to_rank(interior);
    I->add_arc(interior, one);

    // Create set to track nodes.
    std::map<std::vector<int>, Node*> nodes;

    // Main loop.
    std::vector<int> verts(n_verts);
    Eigen::VectorXi color(n_verts);
    for (int i = 0; i < n_verts; i++) {
        verts[i] = i;
    }
    for (int k = d + 1; k > 1; k--) {
        for (Node* H : I->rank(k)) {
            std::vector<int> V_H = difference_sorted(verts, H->_vertices);
            color.setZero();
            for (int v : V_H) {
                std::vector<int> h_v = H->_vertices;
                h_v = insert_sorted(h_v, v);
                std::vector<int> G = closure(h_v, V, F);
                // std::cout << "hv " << H->_vertices << std::endl;
                // std::cout << "in " << h_v << std::endl;
                // std::cout << "cl " << G << std::endl;
                // std::cout << "closure ";
                // print_vector_int(G);
                // std::cout << std::endl;
                if (G.size() == 0) {
                    continue;
                }
                // std::cout << "rank " << k-1 << std::endl;
                // std::cout << "node " << G << std::endl;
                color[v] = 1;
                std::vector<int> W = difference_sorted(G, H->_vertices);
                for (int w : W) {
                    if (w == v) {
                        continue;
                    }
                    if (color[w] >= 0) {
                        color[v] = -1;
                        break;
                    }
                }
                if (color[v] == 1) {
                    Node* u;
                    if (nodes.count(G) == 0) {
                        // Create node.
                        u = I->make_node(k-1);
                        u->_vertices = G;
                        I->add_node_to_rank(u);
                        nodes.insert(std::make_pair(G, u));
                    }
                    else {
                        u = nodes[G];
                    }
                    // Add arcs.
                    I->add_arc(u, H);
                }
            }
        }
    }

    // Add origin.
    Node* origin = I->make_node(0);
    origin->_key = std::string(n, '0');
    I->add_node_to_rank(origin);
    for (Node* u : I->rank(1)) {
        I->add_arc(origin, u);
    }
    origin->_vertices = verts;

    // Add zero node.
    Node* zero = I->make_node(-1);
    zero->_key = std::string(n, 'z');
    zero->sign_vector = "{0}";
    I->add_node_to_rank(zero);
    for (Node* u : I->rank(0)) {
        I->add_arc(zero, u);
    }

    // for (Node* u : I->_nodes) {
    //     std::cout << "(";
    //     for (int i = 0; i < u->_vertices.size(); i++) {
    //         std::cout << u->_vertices[i] << " ";
    //     }
    //     std::cout << ")" <<  std::endl;
    // }

    // Update positions and sign vectors.
    for (Node* u : I->_nodes) {
        if (u->rank == -1 || u->rank >= d + 1) {
            continue;
        }
        u->position.resize(n);
        u->position.setOnes();
        u->sign_vector.resize(n, '-');
        for (int i = 0; i < u->_vertices.size(); i++) {
            int v = u->_vertices[i];
            u->position[v] = 0;
            u->sign_vector[v] = '0';
        }
        // For interior point computation.
        u->_key = u->sign_vector;
        // std::cout << u->sign_vector << std::endl;
    }
    // for (Node* u : I->_nodes) {
    //     std::cout << *u << "\n" << std::endl;
    // }

    // Update interior points.
    for (int k = 0; k <= d; k++) {
        for (Node* u : I->rank(k)) {
            // u->update_interior_point(eps);
        }
    }

    return I;
}

Arc::Arc() 
    : dst(nullptr), _dst_arc(nullptr), 
    _next(nullptr), _prev(nullptr) {}

ArcList::ArcList() : _begin(nullptr), _size(0) {}

void IncidenceGraph::add_arc(Node* src, Node* dst,
                             Arc* arc_src, Arc* arc_dst) {
    // Get corresponding arc-list of source (sub) and destination (super) nodes.
    ArcList* src_list = &src->superfaces;
    ArcList* dst_list = &dst->subfaces;

    // Allocate arcs from memory pool.
    if (!arc_src) {
        arc_src = _arc_pool.malloc();
        _num_arcs_created += 1;
    }
    if (!arc_dst) {
        arc_dst = _arc_pool.malloc();
        _num_arcs_created += 1;
    }

    // Create arcs.
    arc_src->dst = dst;
    arc_src->_dst_arc = arc_dst;
    arc_dst->dst = src;
    arc_dst->_dst_arc = arc_src;

    // Add arcs.
    src_list->_add_arc(arc_src);
    dst_list->_add_arc(arc_dst);
}

void ArcList::_add_arc(Arc* arc) {
    arc->_prev = nullptr;
    if (!_begin) {
        arc->_next = nullptr;
    } else {
        _begin->_prev = arc;
        arc->_next = _begin;
    }
    _begin = arc;
    _size += 1;
}

void IncidenceGraph::remove_arc(Arc* arc) {
    Node* src = arc->_dst_arc->dst;
    Node* dst = arc->dst;

    // Get corresponding arc-lists.
    ArcList* src_list;
    ArcList* dst_list;
    if (src->rank > dst->rank) {
        src_list = &src->subfaces;
        dst_list = &dst->superfaces;
    } else if (src->rank < dst->rank) {
        src_list = &src->superfaces;
        dst_list = &dst->subfaces;
    } else {
        MODUS_ASSERT(false);
    }

    // Remove arcs.
    src_list->_remove_arc(arc);
    dst_list->_remove_arc(arc->_dst_arc);
}

void ArcList::_remove_arc(Arc* arc) {
    _size -= 1;

    // Connect previous and next arcs.
    if (arc->_prev) {
        arc->_prev->_next = arc->_next;
    }
    if (arc->_next) {
        arc->_next->_prev = arc->_prev;
    }

    // Update beginning index.
    if (_begin == arc) {
        _begin = arc->_next;
    }
}

ArcListIterator ArcList::begin() {
    if (!_begin) {
        return end();
    }
    return ArcListIterator(this, _begin);
}

ArcListIterator ArcList::end() {
    return ArcListIterator(this, nullptr);
}

ArcListIterator::ArcListIterator(ArcList* arc_list, Arc* arc) {
    this->arc = arc;
    this->arc_list = arc_list;
}

ArcListIterator& ArcListIterator::operator++() {
    if (!this->arc) {
        throw std::runtime_error("Increment a past-the-end iterator");
    } else if (arc->_next == nullptr) {
        this->arc = nullptr;
    } else {
        this->arc = arc->_next;
    }
    return *this;
}

ArcListIterator ArcListIterator::operator++(int n) {
    ArcListIterator retval = *this; 
    ++(*this); 
    return retval;
}

Node* ArcListIterator::operator*() const {
    return arc->dst;
}

Node** ArcListIterator::operator->() const {
    return &arc->dst;
}

bool operator==(const ArcListIterator& lhs, const ArcListIterator& rhs) {
    return lhs.arc == rhs.arc && lhs.arc_list == rhs.arc_list;
}

bool operator!=(const ArcListIterator& lhs, const ArcListIterator& rhs) {
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& out, const Arc& arc) {
    out << "arc: dst = " << arc.dst << std::endl
        << " dst idx = " << arc._dst_arc << std::endl
        << "    next = " << arc._next << std::endl
        << "    prev = " << arc._prev;
    return out;
}

Node::Node(int k)
    : rank(k), _id(-1), _color(COLOR_AH_WHITE),
      _black_bit(0), _sign_bit_n(0), _sign_bit(0), _graph(nullptr), feasible(false), 
      _side(SIDE_AH_UNSET)

{
    this->_grey_subfaces.clear();
    this->_black_subfaces.clear();
    this->_key.clear();

    this->interior_point.resize(0);
    this->position.resize(0);
    this->sign_vector.clear();
}

void Node::reset() {
    rank = -100;
    _id = -1;
    _color = COLOR_AH_WHITE;
    _black_bit = 0;
    _sign_bit_n = 0;
    _sign_bit = 0;
    _grey_subfaces.clear();
    _black_subfaces.clear();
    _key.clear();
    interior_point.resize(0);
    position.resize(0);
    sign_vector.clear();
}

std::ostream& operator<<(std::ostream& out, Node& node) {
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    out << "node: " << node._id << "\n"
        << "rank: " << (int) node.rank << "\n"
        << " key: " << node._key << "\n"
        << "  sv: " << node.sign_vector << "\n"
        << "  pt:" << std::fixed << std::setprecision(6) << std::setfill(' ') 
        << node.interior_point.transpose().format(CleanFmt) << "\n"
        << "#sub: " << node.subfaces.size() << "\n"
        << " sub: ";
    bool first = true;
    for (Node* f : node.subfaces) {
        if (first) {
            out << f->sign_vector << std::endl;
            first = false;
        } else {
            out << "      " << f->sign_vector << std::endl;
        }
    }
    out 
    << "#sup: " << node.superfaces.size() << "\n"
    << " sup: ";
    first = true;
    for (Node* f : node.superfaces) {
        if (first) {
            out << f->sign_vector;
            first = false;
        } else {
            out << std::endl << "      " << f->sign_vector;
        }
    }
    return out;
}

void Node::update_interior_point(double eps) {
    if (this->rank == -1 || this->rank == this->graph()->dim() + 1) {
        return;
    }

    // Case: Vertex
    else if (this->rank == 0) {
        // Solve linear equations for unique intersection point.
        Eigen::VectorXi idx;
        arg_where(this->_key, '0', idx);
        if (DEBUG) {
            std::cout << this->_key << std::endl;
            std::cout << idx.transpose() << std::endl;
            std::cout << _graph->A << std::endl;
            // MODUS_ASSERT(idx.size() == _graph->A.cols());
        }
        Eigen::MatrixXd A0;
        Eigen::VectorXd b0;
        GetRows(_graph->A, idx, A0);
        GetRows(_graph->b, idx, b0);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
        qr.setThreshold(eps);
        qr.compute(A0);
        this->interior_point = qr.solve(b0);
    }

    // Case: Edge
    else if (this->rank == 1) {

        // Case: 2 vertices: Average interior points of vertices.
        if (this->subfaces.size() == 2) {
            auto iter = this->subfaces.begin();
            Node* v0 = *iter++;
            Node* v1 = *iter;
            this->interior_point = (v0->interior_point + v1->interior_point) / 2.0;
        }

        // Case: 1 vertex: Pick an interior point along the unbounded edge
        // (ray).
        else if (this->subfaces.size() == 1) {
            // Get vertex.
            Node* v = *this->subfaces.begin();
            // Get edge key up to vertex key length.
            std::string e_sv = this->_key.substr(0, v->_key.size());
            if (DEBUG) {
                MODUS_ASSERT(v->_key.size() <= this->_key.size());
            }
            // Compute edge direction.
            Eigen::VectorXi id0;
            arg_where(this->_key, '0', id0);
            Eigen::MatrixXd A0 = GetRows(this->graph()->A, id0);
            Eigen::MatrixXd kernel = 10 * kernel_basis(A0, eps);
            if (kernel.cols() > 1) {
                // Select the column which achieves the lowest residual.
                Eigen::MatrixXd r = (A0 * kernel).transpose();
                double d_min = std::numeric_limits<double>::infinity();
                int i_min = 0;
                for (int i = 0; i < r.rows(); i++) {
                    double d = r.row(i).norm();
                    if (d < d_min) {
                        d_min = d;
                        i_min = i;
                    }
                }
                kernel = kernel.col(i_min);

                if (DEBUG) {
                    std::cout << "Extra kernel bases" << std::endl;
                    std::cout << "Selecting minimum residual basis" << std::endl;
                    std::cout << kernel << std::endl;
                }
            } else if (kernel.cols() == 0) {
                std::cout << "kernel 0" << std::endl;
            }
            if (DEBUG) {
                std::cout << *this << std::endl;
                std::cout << "key " << _key << std::endl;
                std::cout << "id0 " << id0.transpose() << std::endl;
                std::cout << "A0\n" << A0 << std::endl;
                std::cout << "kernl\n" << kernel << std::endl;
                std::cout << "A0*K"<< (A0 * kernel).transpose() << std::endl;
                MODUS_ASSERT(kernel.cols() == 1);
            }
            // Get vertex key.
            std::string v_sv = v->_key;
            if (DEBUG) {
                std::cout << e_sv << std::endl;
                std::cout << v_sv << std::endl;
                std::cout << "A * kernel\n" << (this->graph()->A * kernel).transpose() << std::endl;
            }
            // Find an index where edge and vertex keys disagree.
            Eigen::VectorXi idx;
            arg_not_equal(e_sv, v_sv, idx);
            int i = idx[0];
            if (DEBUG) {
                MODUS_ASSERT(e_sv[i] != '0');
                MODUS_ASSERT(v_sv[i] == '0');
            }
            // Determine ray direction.
            Eigen::MatrixXd A_i = this->graph()->A.row(i);
            Eigen::VectorXd b_i = this->graph()->b.row(i);
            double dot = (A_i * (v->interior_point + kernel) - b_i).sum();
            double sgn = e_sv[i] == '+' ? 1 : -1;
            if (dot * sgn > 0) {
                this->interior_point = v->interior_point + kernel;
            } else {
                this->interior_point = v->interior_point - kernel;
            }
        } else {
            std::cout << "      edge: " << this->_key << std::endl;
            std::cout << "# subfaces: " << this->subfaces.size() << std::endl;
            std::cout << "  subfaces: ";
            for (Node* f : subfaces) {
                std::cout << f->_id << " ";
            }
            std::cout << std::endl;
            MODUS_ASSERT(false);
        }
    }

    // Case: k-face with k >= 2: Average interior points of vertices.
    else {
        this->interior_point.resize(this->graph()->A.cols());
        this->interior_point.setZero();
        for (Node* f : subfaces) {
            this->interior_point += f->interior_point;
        }
        this->interior_point /= this->subfaces.size();
    }
    if (DEBUG) {
        // std::cout << *this << std::endl << std::endl;

        this->update_sign_vector(eps);

        std::cout << "\nTESTING\n" << *this << std::endl << std::endl;

        // Assert superfaces obey partial order.
        for (Node* f : superfaces) {
            if (DEBUG) {
                MODUS_ASSERT(f->rank == this->rank + 1);
            }
        }

        // Assert subfaces obey partial order.
        for (Node* f : subfaces) {
            if (DEBUG) {
                MODUS_ASSERT(f->rank == this->rank - 1);
            }
            if (f->rank == -1 || f->rank == _graph->A.cols()) {
                continue;
            }
            f->update_sign_vector(eps);
            int order = partial_order(f->sign_vector, this->sign_vector);
            if (!(order & LESS_THAN_EQUAL)) {
                // Print sign vectors of subfaces.
                // std::cout << "sv:\n" << this->sign_vector << std::endl;
                // std::cout << "subfaces:" << std::endl;
                // for (int i : this->subfaces) {
                //     NodePtr f = _graph->node(i);
                //     std::cout << f->sign_vector << " " << f->_id << std::endl;
                // }
                MODUS_ASSERT(false);
            }
        }

        // Assert sign vector matches key.
        this->update_sign_vector(eps);
        bool match = this->_key == this->sign_vector.substr(0, this->_key.size());
        if (!match) {
            std::cout << *this->subfaces._begin << std::endl;
            std::cout << *this << std::endl;
        }
        MODUS_ASSERT(match);
    }
}

void Node::update_position(double eps) {
    if (DEBUG) {
        // std::cout << "rank: " << (int) this->rank << std::endl;
        // std::cout << _graph->A << std::endl;
        // std::cout << _graph->b << std::endl;
        // std::cout << "int pt\n" << this->interior_point << std::endl;
    }
    if (this->rank == -1 || this->rank == this->_graph->dim() + 1) {
        return;
    } else {
        Eigen::VectorXd res = _graph->A * this->interior_point - _graph->b;
        get_position(res, this->position, eps);
    }
}

void Node::update_sign_vector(double eps) {
    if (this->rank == -1 || this->rank == this->_graph->dim() + 1) {
        return;
    }
    update_position(eps);
    sign_vector.clear();
    for (int i = 0; i < this->position.size(); i++) {
        if (this->position[i] == 0) {
            sign_vector.push_back('0');
        } else if (this->position[i] == 1) {
            sign_vector.push_back('+');
        } else if (this->position[i] == -1) {
            sign_vector.push_back('-');
        } else {
            MODUS_ASSERT(false);
        }
    }
}

IncidenceGraph::IncidenceGraph(int d) 
    : _num_nodes_created(0), _num_arcs_created(0)
{
    this->_lattice.resize(d + 3);
    for (int i = 0; i < d + 3; i++) {
        this->_lattice[i].clear();
    }
}

IncidenceGraph::~IncidenceGraph() {
}

int IncidenceGraph::dim() {
    return this->_lattice.size()-3;
}

int IncidenceGraph::num_k_faces(int k) {
    return this->rank(k).size();
}

int IncidenceGraph::num_incidences() {
    return -1;
}

int IncidenceGraph::num_ranks() {
    return this->_lattice.size();
}

int IncidenceGraph::num_proper_ranks() {
    return this->_lattice.size() - 2;
}

int IncidenceGraph::rank() {
    return this->_lattice.size() - 2;
}

int IncidenceGraph::rank_facet() {
    return this->_lattice.size() - 3;
}

int IncidenceGraph::max_zero_rank() {
    for (int i = 0; i < _lattice.size(); i++) {
        if (_lattice[i].size() > 1) {
            return i - 2;
        }
    }
    return -std::numeric_limits<int>::infinity();
}

int IncidenceGraph::min_one_rank() {
    for (int i = _lattice.size() - 1; i >= 0; i--) {
        if (_lattice[i].size() > 1) {
            return i;
        }
    }
    return std::numeric_limits<int>::infinity();
}

void IncidenceGraph::set_hyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    this->A = A;
    this->b = b;
}

void IncidenceGraph::add_hyperplane(const Eigen::VectorXd& a, double d) {
    int n = this->A.rows();
    int m = this->A.cols();
    Eigen::MatrixXd A(n + 1, m);
    A << this->A, a.transpose();
    Eigen::VectorXd b(n + 1);
    b << this->b, d;
    this->A = A;
    this->b = b;
    // if (DEBUG) {
    //     std::cout << "a\n" << a << std::endl;
    //     std::cout << "A\n" << this->A << std::endl;
    //     std::cout << "b\n" << d << std::endl;
    //     std::cout << "b\n" << this->b << std::endl;
    // }
}

void IncidenceGraph::update_positions(double eps) {
    for (int i = 0; i < this->_nodes.size(); i++) {
        this->_nodes[i]->update_position(eps);
    }
}

void IncidenceGraph::update_sign_vectors(double eps) {
    for (int i = 0; i < this->_nodes.size(); i++) {
        this->_nodes[i]->update_sign_vector(eps);
    }
}

Positions IncidenceGraph::get_positions() {
    Positions P;
    for (int i = 0; i < this->_nodes.size(); i++) {
        if (this->_nodes[i]->position.size() != 0) {
            P.push_back(this->_nodes[i]->position);
        }
    }
    return P;
}

SignVectors IncidenceGraph::get_sign_vectors() {
    SignVectors S;
    for (int i = 0; i < this->_nodes.size(); i++) {
        if (!this->_nodes[i]->sign_vector.empty()) {
            S.push_back(this->_nodes[i]->sign_vector);
        }
    }
    return S;
}

SignVectors IncidenceGraph::get_proper_sign_vectors() {
    SignVectors S;
    for (int k = 0; k < this->num_proper_ranks(); k++) {
        for (Node* u : this->rank(k)) {
            S.push_back(u->sign_vector);
        }
    }
    return S;
}

SignVectors IncidenceGraph::get_sign_vectors(const std::string& cs_mode) {
    SignVectors S;
    for (int i = 0; i < this->_nodes.size(); i++) {
        const std::string& sv = this->_nodes[i]->sign_vector;
        if (sv.substr(0, cs_mode.size()) == cs_mode) {
            S.push_back(sv);
        }
    }
    return S;
}

inline bool
is_aligned(const void * ptr, std::uintptr_t alignment) noexcept {
    auto iptr = reinterpret_cast<std::uintptr_t>(ptr);
    return !(iptr % alignment);
}

Node* IncidenceGraph::make_node(int k) {
    // Node* node = _node_pool.malloc();
    Node* node = _node_pool.construct(k);
    // node->reset();
    node->rank = k;
    node->_id = this->_nodes.size();
    node->_graph = this;
    this->_nodes.push_back(node);
    if (DEBUG) {
        // Assert cacheline aligned.
    }
    return node;
}

void IncidenceGraph::add_node_to_rank(Node* node) {
    this->rank(node->rank).push_back(node);
}

void IncidenceGraph::remove_node_from_rank(Node* node) {
    // Remove node.
    auto r = rank(node->rank);
    auto p = std::find(r.begin(), r.end(), node);
    r.erase(p);
}

Rank& IncidenceGraph::rank(int k) {
    return this->_lattice[k+1];
}