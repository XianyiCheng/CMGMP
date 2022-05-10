#include <modus/modes/geometry/convex_hull.hpp>
#include <modus/modes/geometry/incidence_graph.hpp>
// #include <glog/logging.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>

#include <libqhullcpp/QhullLinkedList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullPoints.h>
#include <libqhullcpp/QhullHyperplane.h>

using orgQhull::Qhull;
using orgQhull::QhullPoint;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullHyperplane;

static int DEBUG=0;
static int PROFILE=0;


Eigen::MatrixXi convex_hull_1d(const Eigen::MatrixXd& pts, double eps) {
    // Take convex hull, i.e. min + max points, in 1D.
    int i_min = pts.minCoeff();
    int i_max = pts.maxCoeff();

    // Find min and max points.
    double pt_min = pts.minCoeff();
    double pt_max = pts.maxCoeff();
    std::vector<int> min_idx;
    std::vector<int> max_idx;
    for (int i = 0; i < pts.cols(); i++) {
        if (std::abs(pt_min - pts(0,i)) < eps) {
            min_idx.push_back(i);
        }
        if (std::abs(pt_max - pts(0,i)) < eps) {
            max_idx.push_back(i);
        }
    }

    // Create vertex-facet matrix. Find duplicate points.
    Eigen::MatrixXi M(pts.cols(), 2);
    M.setZero();
    for (int i = 0; i < min_idx.size(); i++) {
        M(min_idx[i], 0) = 1;
    }
    for (int i = 0; i < max_idx.size(); i++) {
        M(max_idx[i], 1) = 1;
    }

    // std::cout << "M\n" << M << std::endl;

    return M;
}

/**
 * @brief Computes the convex hull of the input points and returns the vertex
 * facet incidence matrix. 
 *
 * Input has size d×n, where d is the dimension and n is the number of points.
 *
 * @param pts 
 * @return Eigen::MatrixXi 
 */
Eigen::MatrixXi convex_hull(const Eigen::MatrixXd& pts, double eps) {
    // Get dimensions.
    int n = pts.cols();
    int d = pts.rows();

    // Handle d=1 case.
    if (d == 1) {
        return convex_hull_1d(pts, eps);
    }

    // Copy to coordinate type.
    std::vector<coordT> data;
    for (int i = 0; i < pts.size(); i++) {
        data.push_back(pts.data()[i]);
    }

    // Take convex hull.
    Qhull q;
    q.setFactorEpsilon(1e-7);
    std::stringstream buffer;
    buffer << "Fv C" << std::fixed << eps;
    q.runQhull("convex", d, n, data.data(), buffer.str().c_str());
    // q.checkIfQhullInitialized();

    // Build vertex-facet incidence matrix. 
    std::ostringstream oss;
    q.setOutputStream(&oss);
    q.outputQhull();
    Eigen::MatrixXi M;
    std::istringstream qss(oss.str());
    std::string line, token;
    std::getline(qss, line);
    int n_f = std::stoi(line);
    M.setZero(n, n_f);
    for (int i = 0; i < n_f; i++) {
        std::getline(qss, line);
        // std::cout << line << std::endl;
        std::istringstream fv(line);
        std::getline(fv, token, ' ');
        while(std::getline(fv, token, ' ')) {
            int j = std::stoi(token);
            M(j, i) = 1;
        }
    }

    // 
    QhullPoints points = q.points();
    QhullPointsIterator iter(points);
    while(iter.hasNext()) {
        QhullPoint pt = iter.next();
        std::vector<double> ptv = pt.toStdVector();
        for (int i = 0; i < ptv.size(); i++) {
            // std::cout << ptv[i] << " " << std::endl;
        }
        std::string qstr(1000, 0);
        // pt.print(qstr.data());
        // std::cout << "pt" << qstr << std::endl;
    }

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "\t\t[", "]");
    // DLOG(INFO) << "\n\t" << "convex hull vf-matrix"
    //            << "\n" << M.format(CleanFmt)
    //            << "\n\t" << "pts aff"
    //            << std::fixed << std::setprecision(6) << std::setfill(' ')
    //            << "\n" << pts.format(CleanFmt);

    Eigen::MatrixXi M_;
    M_.setZero(n, n_f);
    q.defineVertexNeighborFacets();
    QhullFacetList F(q.beginFacet(), q.endFacet());
    int i_f = 0;
    for (QhullFacetListIterator iter = F; iter.hasNext(); ) {
        QhullFacet f = iter.next();
        QhullHyperplane h = f.hyperplane();
        std::vector<double> hvec = h.toStdVector();
        Eigen::VectorXd h_(hvec.size()-1);
        for (int i = 0; i < hvec.size()-1; i++) {
            h_[i] = hvec[i];
        }
        // std::cout << "h: " << h_.transpose().format(CleanFmt) << " | " << 
        // h.offset() << std::endl;
        // std::cout << std::fixed << std::setprecision(6) << std::setfill(' ') << 
        //     "h" << ((h_.transpose() * pts).array() + h.offset()).format(CleanFmt) << std::endl;
        // std::cout << std::fixed << std::setprecision(6) << std::setfill(' ') << 
        //     "h" << get_sign_vector(Eigen::VectorXd((h_.transpose() * pts).array() + h.offset()), eps)
        //     << std::endl;
        std::string sv = get_sign_vector(Eigen::VectorXd((h_.transpose() * pts).array() + h.offset()), eps);
        for (int i = 0; i < sv.size(); i++) {
            if (sv[i] == '0') {
                M_(i, i_f) = 1;
            }
        }
        i_f ++;
    }
    // std::cout << "M_\n" << M_ << std::endl;

    return M_;
}