#pragma once
#include <memory>
#include <modus/collision/contact.hpp>
#include <modus/kinematics/jacobian.hpp>
#include <modus/modes/geometry/incidence_graph.hpp>


enum {
    MODES_INTERIOR_POINT = 1,
    MODES_DEBUG = 32,
    MODES_PROFILE = 64,
};

struct CSModesOptions {
    int flags;
};

typedef std::shared_ptr<CSModesOptions> CSModesOptionsPtr;

struct SSModesOptions {

};

typedef std::shared_ptr<SSModesOptions> SSModesOptionsPtr;

struct ModeEnumerationOptions {
    int             flags = 0;
    Eigen::VectorXd interior_point;
    int             PROFILE = 0;
};

typedef std::shared_ptr<ModeEnumerationOptions> ModeEnumerationOptionsPtr;

/**
 * @brief Enumerate contacting/separating modes. Normal velocities constraints
 * should follow the convention:
 *      Nx - d <= 0.
 *
 * @param N 
 * @param d 
 * @param eps 
 * @param options 
 * @return IncidenceGraph* 
 */
IncidenceGraph* enumerate_cs_modes(const Eigen::MatrixXd& N, 
                                   const Eigen::VectorXd& d, double eps,
                                   ModeEnumerationOptions* options=nullptr);

/**
 * @brief Enumerate sliding/sticking modes given a contacting/separating mode.
 * 
 * @param N 
 * @param d 
 * @param T 
 * @param cs_mode
 * @param eps 
 * @param options 
 * @return IncidenceGraph* 
 */
IncidenceGraph* enumerate_ss_modes(const Eigen::MatrixXd& N,
                                   const Eigen::VectorXd& d,
                                   const Eigen::MatrixXd& T, 
                                   std::string cs_mode, double eps,
                                   ModeEnumerationOptions* options=nullptr);

/**
 * @brief Enumerate sliding/sticking modes given a contacting/separating mode
 * graph.
 *
 * @param cs_graph 
 * @param N 
 * @param d 
 * @param T 
 * @param eps 
 * @param options 
 * @return std::vector<IncidenceGraph*> 
 */
std::vector<IncidenceGraph*> 
enumerate_ss_modes(IncidenceGraph* cs_graph, const Eigen::MatrixXd& N, 
                   const Eigen::VectorXd& d, const Eigen::MatrixXd& T, double eps,
                   ModeEnumerationOptions* options=nullptr);


namespace modus
{

IncidenceGraph* EnumerateModes(const Eigen::MatrixXd& N, const Eigen::MatrixXd& T,
                               double eps, ModeEnumerationOptions* options);

IncidenceGraph* EnumerateCSModes(const Eigen::MatrixXd& N, double eps, 
                                 ModeEnumerationOptions* options);

IncidenceGraph* EnumerateSSModes(const Eigen::MatrixXd& N, const Eigen::MatrixXd& T,
                                 const std::string& cs_mode, double eps, 
                                 ModeEnumerationOptions* options);

}