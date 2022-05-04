#pragma once
#include <modus/common/eigen.hpp>
#include <memory>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include <boost/pool/object_pool.hpp>


class Node;
typedef std::shared_ptr<Node> NodePtr;

class IncidenceGraph;
typedef std::shared_ptr<IncidenceGraph> IncidenceGraphPtr;

/**
 * @brief Colors for hyperplane arrangement algorithm.
 * f = face (relative interior)
 * h = hyperplane
 * cl = closure
 */
enum {
    COLOR_AH_WHITE   = 0,   /// if cl(f) ∩ h = ∅
    COLOR_AH_PINK    = 1,   /// if cl(f) ∩ h ≠ ∅, f ∩ h = ∅
    COLOR_AH_RED     = 2,   /// if f ∩ h ≠ ∅, f ⊈ h
    COLOR_AH_CRIMSON = 3,   /// if f ⊆ h
    COLOR_AH_GREY    = 4,   /// if cl(f) ∩ h ≠ ∅, f ∩ h = ∅
    COLOR_AH_BLACK   = 5,   /// if f ⊆ h
    COLOR_AH_GREEN   = 6,   /// if some subface is non-white
};

enum {
    SIDE_AH_UNSET   = 0,    /// if f ⊆ ?
    SIDE_AH_ABOVE   = 1,    /// if f ⊆ h₊
    SIDE_AH_BELOW   = 2,    /// if f ⊆ h₋
    SIDE_AH_RED     = 3,    /// if f ∩ h ≠ ∅, f ⊈ h
    SIDE_AH_CRIMSON = 4,    /// if f ⊆ h
};

std::string get_color_ah_string(int color);

typedef Eigen::VectorXi Position;
typedef std::vector<Position> Positions;
typedef std::string SignVector;
typedef std::vector<SignVector> SignVectors;

/**
 * @brief Get the position p(x) ∈ [-1, 0, 1], with p(x) = 0 iff x ∈ [-ε, ε].
 * 
 * @param x 
 * @param eps 
 * @return int 
 */
int  get_position(double x, double eps);

/**
 * @brief Get the sign s(x) ∈ [-, 0, +], with s(x) = 0 iff x ∈ [-ε, ε].
 * 
 * @param x 
 * @param eps 
 * @return char 
 */
char get_sign(double x, double eps);

/**
 * @brief Get the position vector p(v) ∈ [-1, 0, 1]ⁿ.
 * 
 * @param v 
 * @param pos 
 * @param eps 
 */
void get_position(const Eigen::VectorXd& v, Eigen::VectorXi& pos, double eps);

/**
 * @brief Get the position vector p(v) ∈ [-1, 0, 1]ⁿ.
 * 
 * @param v 
 * @param eps 
 * @return Eigen::VectorXi 
 */
Eigen::VectorXi get_position(const Eigen::VectorXd& v, double eps);

/**
 * @brief Get the sign vector s(v) ∈ [-, 0, +]ⁿ.
 * 
 * @param v 
 * @param sv 
 * @param eps 
 */
void get_sign_vector(const Eigen::VectorXd& v, std::string& sv, double eps);

/**
 * @brief Get the sign vector s(v) ∈ [-, 0, +]ⁿ.
 * 
 * @param v 
 * @param eps 
 * @return std::string 
 */
std::string get_sign_vector(const Eigen::VectorXd& v, double eps);

/**
 * @brief Return indices where s(v) is equal to sign z.
 * 
 * @param sv 
 * @param s 
 * @param idx 
 */
void arg_where(const std::string& sv, char z, Eigen::VectorXi& idx);

/**
 * @brief Return indices where x is equal 0.
 * 
 * @param x 
 * @param eps 
 * @param idx 
 */
void arg_where_zero(const Eigen::VectorXd& x, double eps, Eigen::VectorXi& idx);

/**
 * @brief Return indices where x is not equal to 0.
 * 
 * @param x 
 * @param eps 
 * @param idx 
 */
void arg_where_nonzero(const Eigen::VectorXd& x, double eps, Eigen::VectorXi& idx);

/**
 * @brief Return indices where a == b. Must be of same length.
 * 
 * @param a 
 * @param b 
 * @param idx 
 */
void arg_equal(const std::string& a, const std::string& b, Eigen::VectorXi& idx);

/**
 * @brief Return indices where a != b. Must be of same length.
 * 
 * @param a 
 * @param b 
 * @param idx 
 */
void arg_not_equal(const std::string& a, const std::string& b, Eigen::VectorXi& idx);
// bool less_than(const Eigen::VectorXi& a, const Eigen::VectorXi& b);

/**
 * @brief Returns a ∩ b. Inputs must be in sorted order.
 * 
 * @param a 
 * @param b 
 * @return std::vector<int>
 */
std::vector<int> intersect_sorted(const std::vector<int>& a, 
                                  const std::vector<int>& b);

/**
 * @brief Returns c = a / b. Inputs must be in sorted order.
 * 
 * @param a 
 * @param b 
 * @return std::vector<int>
 */
std::vector<int> difference_sorted(const std::vector<int>& a, 
                                   const std::vector<int>& b);

/**
 * @brief Insert v into h. h must be in sorted order.
 * 
 * @param h 
 * @param v 
 * @return std::vector<int> 
 */
std::vector<int> insert_sorted(const std::vector<int>& h, int v);

/**
 * @brief Returns ∩ B(a) over a ∈ A. Inputs must be in sorted order.
 * 
 * @param A 
 * @param B 
 * @return std::vector<int> 
 */
std::vector<int> intersect_preimages(const std::vector<int>& A, 
                                     const std::vector<std::vector<int> >& B);

/**
 * @brief Return the closure S' = F(V(S)) where A(B) denotes preimage
 * intersection. Inputs must be in sorted order.
 *
 * @param S 
 * @param V 
 * @param F 
 * @return std::vector<int> 
 */
std::vector<int> closure(const std::vector<int>& S, 
                         const std::vector<std::vector<int> >& V,
                         const std::vector<std::vector<int> >& F);

/**
 * @brief Constructs the incidence graph (lattice) of an h-polyhedral pointed
 * cone from its vertex figure's (h-cone/origin) vertex-facet incidence matrix.
 *
 * The input dimension d is of the vertex figure. The h-cone will have dimension
 * d + 1.
 *
 * @param M     vertex-facet incidence matrix
 * @param n     number of points
 * @param d     dimension
 * @param eps   tolerance
 * @return IncidenceGraph* 
 */
IncidenceGraph* build_graph_cone(const Eigen::MatrixXi& M, int n, int d, double eps);

/**
 * @brief Enum for partial order.
 */
enum {
    STRICTLY_LESS       = 1,
    STRICTLY_GREATER    = 2,
    EQUAL               = 4,
    LESS_THAN_EQUAL     = 5,
    GREATER_THAN_EQUAL  = 6,
    INCOMPARABLE        = 8
};

/**
 * @brief Return the partial order of lhs ? rhs.
 * 
 * @param lhs sign vector
 * @param rhs sign vector
 * @return int 
 */
int partial_order(const std::string& lhs, const std::string& rhs);

/**
 * @brief An arc (directed edge) between two faces of the incidence graph.
 * 
 */
struct Arc {
    Node* dst;          /// Node at the end of this arc.
    Arc*  _dst_arc;
    Arc*  _next;
    Arc*  _prev;

    Arc();

    friend std::ostream& operator<<(std::ostream& out, const Arc& arc);
};

class ArcListIterator;

typedef std::vector<Arc> Arcs;

/**
 * @brief Doubly linked link of arcs for a node, i.e. subfaces or superfaces.
 * 
 */
class ArcList {
public:
    Arc* _begin;
    int  _size;

    ArcList();

    int size() { return _size; }
    ArcListIterator begin();
    ArcListIterator end();

    friend IncidenceGraph;

    void _add_arc(Arc* arc);
    void _remove_arc(Arc* arc);
};

/**
 * @brief Forward iterator for an arc list.
 * 
 */
class ArcListIterator {
public:
    using value_type = int;
    using reference = int;
    using iterator_category = std::input_iterator_tag;
    using pointer = int*;
    using difference_type = void;

    Arc*     arc;
    ArcList* arc_list;
    
    class postinc_return {
    public:
        int value;
        postinc_return(int value_) { value = value_; }
        int operator*() { return value; }
    };

    ArcListIterator(ArcList* arc_list, Arc* arc);

    ArcListIterator& operator++();
    ArcListIterator  operator++(int);
    Node*  operator*()  const;
    Node** operator->() const;

    friend bool operator==(const ArcListIterator& lhs, const ArcListIterator& rhs);
    friend bool operator!=(const ArcListIterator& lhs, const ArcListIterator& rhs);
};

/**
 * @brief Class representing a node (face) of the incidence graph (lattice).
 * 
 */
class Node {
public:
    int8_t              rank;               /// Rank of this face.
    int8_t              _color;
    int8_t              _black_bit;
    int8_t              _sign_bit;          
    int8_t              _side;              /// Side of incoming hyperplane.
    int                 _sign_bit_n;
    ArcList             superfaces;         /// Superface arc list.
    ArcList             subfaces;           /// Subface arc list.
    IncidenceGraph*     _graph;
    std::vector<Node*>  _grey_subfaces;
    std::vector<Node*>  _black_subfaces;
    std::string         _key;               // target sign vector
    std::vector<int>    _vertices;          

    int                 _id;
    Eigen::VectorXd     interior_point;     /// A point in its relative interior.
    Position            position;           /// Face position vector.
    SignVector          sign_vector;        /// Face sign vector.
    bool                feasible;           /// Represents a feasible velocity.

    Node(int k);
    void reset();

    IncidenceGraph* graph() { return _graph; }

    void update_interior_point(double eps);
    void update_position(double eps);
    void update_sign_vector(double eps);

    friend std::ostream& operator<<(std::ostream& out, Node& node);
};

typedef std::vector<Node*> Rank;

struct aligned_allocator_mm_malloc_free
{
  typedef std::size_t size_type; //!< An unsigned integral type that can represent the size of the largest object to be allocated.
  typedef std::ptrdiff_t difference_type; //!< A signed integral type that can represent the difference of any two pointers.

  static char * malloc BOOST_PREVENT_MACRO_SUBSTITUTION(const size_type bytes)
  { return static_cast<char *>((_mm_malloc)(bytes, 64)); }
  static void free BOOST_PREVENT_MACRO_SUBSTITUTION(char * const block)
  { (_mm_free)(block); }
};

typedef boost::object_pool<Arc, aligned_allocator_mm_malloc_free> ArcPool;
typedef boost::object_pool<Node, aligned_allocator_mm_malloc_free> NodePool;

/**
 * @brief The incidence graph (lattice) associated with an convex polyhedra or
 * hyperplane arrangement.
 * 
 * Positions (and sign vectors) follow the convention #FIXME
 *  -1 iff Ax - b < 0
 *   0 iff Ax - b = 0
 *  +1 iff Ax - b > 0
 *
 */
class IncidenceGraph : public std::enable_shared_from_this<IncidenceGraph> {
public:
    Eigen::MatrixXd A;                      /// Hyperplane normals.
    Eigen::VectorXd b;                      /// Hyperplane offsets.
    std::vector<Node*> _nodes;
    std::vector<Rank>  _lattice;
    int                _num_nodes_created;
    int                _num_arcs_created;
    NodePool           _node_pool;
    ArcPool            _arc_pool;

    IncidenceGraph(int d);
    ~IncidenceGraph();

    void set_hyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    void add_hyperplane (const Eigen::VectorXd& a, double b);

    int dim();
    int num_k_faces(int k);
    int num_incidences();
    int num_ranks();
    int num_proper_ranks();
    int rank();
    int rank_facet();
    int max_zero_rank();
    int min_one_rank();

    void update_positions(double eps);
    Positions get_positions();

    void update_sign_vectors(double eps);
    SignVectors get_sign_vectors();
    SignVectors get_proper_sign_vectors();
    SignVectors get_sign_vectors(const std::string& cs_mode);

    Node*       node(int id) { return _nodes[id]; }
    Node*  make_node(int k);
    void    add_node_to_rank(Node* node);
    void remove_node_from_rank(Node* node);

    /**
     * @brief Add arc between sub-face and super-face.
     * 
     * @param sub 
     * @param super 
     * @param arc1 
     * @param arc2 
     */
    void add_arc(Node* sub, Node* super,
                 Arc* arc1=nullptr, Arc* arc2=nullptr);   
    // O(1) remove
    void remove_arc(Arc* arc);

    Rank& rank(int k);
};