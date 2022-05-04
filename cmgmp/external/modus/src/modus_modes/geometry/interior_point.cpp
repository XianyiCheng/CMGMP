#include <modus/modes/geometry/interior_point.hpp>
#include <modus/common/linear_algebra.hpp>
#include <modus/common/assert.hpp>
#include <ortools/linear_solver/linear_solver.h>

static int DEBUG=0;


Eigen::VectorXd interior_point(const Eigen::MatrixXd& A,
                               const Eigen::VectorXd& b, double eps) {
                            //    const Eigen::MatrixXd* Aeq,
                            //    const Eigen::VectorXd* beq) {
    using namespace operations_research;

    // Find dimension of halfspaces.
    int n = A.rows();
    int d = A.cols();

    // TODO Reduce A to minimal row rank. (OR-Tools might handle this
    // situation.)

    // Create linear solver.
    MPSolver solver("interior_point", MPSolver::CLP_LINEAR_PROGRAMMING);
    // MPSolver* solver = operations_research::MPSolver::CreateSolver();

    // Create variable x, the interior point.
    std::vector<MPVariable*> x; x.clear();
    double box = 1; // ‖x‖_∞ constraint
    solver.MakeNumVarArray(d, -box, box, "x", &x);

    // Create variable c, the minimum distance from the hyperplanes.
    double inf = std::numeric_limits<double>::infinity();
    MPVariable* c = solver.MakeNumVar(0, inf, "c");

    // DLOG(INFO) << "Number of variables " << solver.NumVariables();

    // Create hyperplane constraints.
    for (int i = 0; i < n; i++) {
        MPConstraint* const ct = solver.MakeRowConstraint(-inf, b[i], std::to_string(i));
        for (int j = 0; j < d; j++) {
            ct->SetCoefficient(x[j], A(i,j));
        }
        ct->SetCoefficient(c, 1);
    }

    // DLOG(INFO) << "Number of constraints " << solver.NumConstraints();

    // Create objective that maximizes the minimum signed distance from the
    // hyperplanes.
    MPObjective* const cost = solver.MutableObjective();
    for (int i = 0; i < d; i++) {
        cost->SetCoefficient(x[i], 0);
    }
    cost->SetCoefficient(c, 1);
    cost->SetMaximization();

    // Solve linear program for strictly interior point.
    solver.Solve();

    // DLOG(INFO) << "Solution:" << std::endl;
    // DLOG(INFO) << "Min dist: " << cost->Value() << std::endl;

    if (DEBUG) {
        std::string model;
        solver.ExportModelAsLpFormat(false, &model);
        std::cout << model << std::endl;
    }

    Eigen::VectorXd x_(d);
    for (int i = 0; i < d; i++) {
        x_[i] = x[i]->solution_value();
    }
    return x_;
}

Eigen::VectorXd modus::InteriorPoint(const Eigen::MatrixXd& N,
                                     const Eigen::MatrixXd& T,
                                     const std::string& cs_mode,
                                     const std::string& ss_mode,
                                     double eps)
{
    // Get various dimensions.
    size_t dimension = N.cols();
    size_t n_contacts = N.rows();
    size_t n_contacting = std::count(cs_mode.begin(), cs_mode.end(), '0');
    MODUS_ASSERT(ss_mode.size() % n_contacting == 0);
    size_t n_sticking = std::count(ss_mode.begin(), ss_mode.end(), '0');
    size_t n_tangent = ss_mode.size() / n_contacting;
    size_t n_total = cs_mode.size() + ss_mode.size();
    
    if (n_contacting == n_contacts && n_sticking == ss_mode.size()) {
        return Eigen::VectorXd::Zero(dimension);
    }

    // Partition hyperplanes into equality constraints and strict inequality
    // constraints.
    Eigen::MatrixXd A(n_total - n_contacting - n_sticking, dimension);
    Eigen::MatrixXd Aeq(n_contacting + n_sticking, dimension);
    size_t ia = 0;
    size_t iae = 0;
    for (size_t i = 0; i < n_contacts; i++) {
        MODUS_ASSERT(cs_mode[i] == '0' || cs_mode[i] == '-');
        if (cs_mode[i] == '0') {
            Aeq.row(iae++) = N.row(i);
        } else {
            A.row(ia++) = N.row(i);
        }
    }
    size_t iss = 0;
    for (size_t i = 0; i < n_contacts; i++) {
        if (cs_mode[i] == '0') {
            for (size_t j = 0; j < n_tangent; j++) {
                switch (ss_mode[iss++]) {
                    case '0':
                        Aeq.row(iae++) = T.row(i*n_tangent + j);
                        break;
                    case '-':
                        A.row(ia++) = T.row(i*n_tangent + j);
                        break;
                    case '+':
                        A.row(ia++) = -T.row(i*n_tangent + j);
                        break;
                    default:
                        MODUS_ASSERT(false);
                        break;
                }
            }
        }
    }
    MODUS_ASSERT(ia == n_total - n_contacting - n_sticking);
    MODUS_ASSERT(iae == n_contacting + n_sticking);
    MODUS_ASSERT(iss == ss_mode.size());

    // Project into nullspace of equality constraints.
    Eigen::MatrixXd K;
    K.setIdentity(dimension, dimension);
    if (n_contacting > 0) {
        K = kernel_basis(Aeq, eps);
        A = A * K;
        MODUS_ASSERT(K.rows() > 0);
        if (K.cols() == 0) {
            std::cout << cs_mode << " " << ss_mode << std::endl;
            std::cout << std::fixed << std::setprecision(5) << Aeq << std::endl;
        }
        MODUS_ASSERT(K.cols() > 0);
    }
    Eigen::VectorXd z;
    z.setZero(A.rows());

    Eigen::VectorXd x = ::interior_point(A, z, eps);

    return K * x;
}