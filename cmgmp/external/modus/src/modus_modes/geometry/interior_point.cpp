#include <iostream>
#include <modus/modes/geometry/interior_point.hpp>
#include <modus/common/linear_algebra.hpp>
#include <modus/common/assert.hpp>

// use lp instead of ortools
#include <glpk.h>

bool lp(const Eigen::VectorXd &C, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu, Eigen::VectorXd *xs, double *optimal_cost) {
  /* declare variables */
  if (xs->rows() <= 0) {
    std::cerr << "[lp] Error: xs is not initialized!!" << std::endl;
    exit(-1);
  }
  glp_prob *lp;
  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.presolve = GLP_OFF;
  parm.msg_lev = GLP_MSG_ERR; // error and warning only
  int *ia, *ja;
  double *ar;
  int rows = A.rows();
  int cols = A.cols();
  int rows_e = Ae.rows();
  int cols_e = Ae.cols();
  // assert(cols_e == cols);

  Eigen::VectorXd xu_expand = Eigen::VectorXd(cols) * nan("");
  Eigen::VectorXd xl_expand = Eigen::VectorXd(cols) * nan("");
  if (xu.rows() > 0) {
    xu_expand = xu;
  }
  if (xl.rows() > 0) {
    xl_expand = xl;
  }

  int size = rows * cols + rows_e * cols_e;
  ia = new int[size + 1000];
  ja = new int[size + 1000];
  ar = new double[size + 1000];

  /**
   * Create problem
   */
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN); // minimization, not maximization

  /**
   * Fill problem
   */
  /* sign and right-hand-side of constraints */
  glp_add_rows(lp, rows + rows_e);
  for (int r = 1; r <= rows; ++r) {
    glp_set_row_bnds(lp, r, GLP_UP, 0.0, b(r-1)); // upper bound, <=
  }
  for (int r = 1; r <= rows_e; ++r) {
    glp_set_row_bnds(lp, rows+r, GLP_FX, be(r-1), be(r-1)); // equality =
  }

  /* cost function */
  glp_add_cols(lp, cols);
  /* variable bounds */
  for (int c = 1; c <= cols; ++c) {
    glp_set_obj_coef(lp, c, C(c-1)); // cost function
    if (std::isfinite(xl_expand(c-1))) {
      if (std::isfinite(xu_expand(c-1))) {
        glp_set_col_bnds(lp, c, GLP_DB, xl_expand(c-1), xu_expand(c-1)); // double bounded
      } else {
        glp_set_col_bnds(lp, c, GLP_LO, xl_expand(c-1), 0.0); // lower-bounded
      }
    } else {
      if (std::isfinite(xu_expand(c-1))) {
        glp_set_col_bnds(lp, c, GLP_UP, 0.0, xu_expand(c-1)); // upper-bounded
      } else {
        glp_set_col_bnds(lp, c, GLP_FR, 0.0, 0.0); // no boundary
      }
    }
  }
  /* fill in coefficient matrix */
  int id = 0;
  for (int r = 1; r <= rows; ++r) {
    for (int c = 1; c <= cols; ++c) {
      id = (r-1)*cols + c;
      ia[id] = r, ja[id] = c, ar[id] = A(r-1, c-1);
    }
  }
  for (int r = 1; r <= rows_e; ++r) {
    for (int c = 1; c <= cols_e; ++c) {
      id = (r + rows - 1)*cols + c;
      ia[id] = r + rows, ja[id] = c, ar[id] = Ae(r-1, c-1);
    }
  }
  glp_load_matrix(lp, id, ia, ja, ar);
  /**
   * solve problem
   */
  // glp_write_prob(lp, 0, "problem.txt");
  // parm.msg_lev = GLP_MSG_ALL;
  glp_simplex(lp, &parm);
  int result = glp_get_status(lp);

  // glp_iptcp iparm;
  // glp_init_iptcp(&iparm);
 
  // // iparm.msg_lev = GLP_MSG_ALL; // error and warning only
  // glp_interior(lp, &iparm);
  // int result = glp_ipt_status(lp);

  /* housekeeping */
  glp_delete_prob(lp);
  glp_free_env();
  delete [] ia;
  delete [] ja;
  delete [] ar;

  if ((result == GLP_OPT) || (result == GLP_FEAS)) {
    // feasible
    *optimal_cost = glp_get_obj_val(lp);
    for (int d = 0; d < cols; ++d) {
      (*xs)(d) = glp_get_col_prim(lp, d + 1);
    }
    // std::cout << "z: " << z << std::endl;
    // std::cout << "solution: " << xs.transpose() << std::endl;
    return true;
  } else {
    return false;
    // std::cout << "Infeasible." << std::endl;
  }

}

static int DEBUG=0;

Eigen::VectorXd interior_point(const Eigen::MatrixXd& A,
                               const Eigen::VectorXd& b, double eps) {

    // Find dimension of halfspaces.
    int n = A.rows();
    int d = A.cols();


    double box = 1; // ‖x‖_∞ constraint

    // TODO Reduce A to minimal row rank.
    Eigen::VectorXd C = Eigen::VectorXd::Zero(d+1);
    C(d) = -1;

    Eigen::MatrixXd A_(n, d+1);
    A_.block(0,0,n,d) = A;
    A_.block(0,d,n,1) = Eigen::VectorXd::Constant(n, 1);

    Eigen::MatrixXd Ae;
    Eigen::VectorXd be;

    Eigen::VectorXd xl = Eigen::VectorXd::Constant(d+1, -box);
    xl(d) = 0;
    Eigen::VectorXd xu = Eigen::VectorXd::Constant(d+1, box);
    xu(d) = 1e+10;

    Eigen::VectorXd xs(d+1); 
    double optimal_cost;
    bool result;
    result = lp(C, A_, b, Ae, be, xl, xu, &xs, &optimal_cost);

    Eigen::VectorXd x_(d);
    for (int i = 0; i < d; i++) {
        x_[i] = xs[i];
    }
    return x_;
}

/*
#include <ortools/linear_solver/linear_solver.h>

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
    std::cout << "From original " << x_.transpose() << std::endl;

    Eigen::VectorXd x1 = interior_point_lp(A, b, eps);
    std::cout << "From lp " << x1.transpose() << std::endl; 
    return x_;
}
*/

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