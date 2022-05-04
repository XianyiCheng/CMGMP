#include <modus/qpcc/solver.hpp>

namespace modus
{

class FullModeSolver : public QPCCSolver {
 public:
  size_t n_contacts_;
  std::vector<std::string> modes_;

  void SetModes(const std::vector<std::string>& modes, size_t n);

  void Solve(QPCCPtr qpcc);

  Eigen::VectorXd GetSolution();
};

}