#pragma once
#include <qpcc/SolverInterface.hpp>


namespace qpcc
{

class AdjacentModesSolver : public QPCCSolver {
public:
    QPCCResultPtr Solve(QPCCProblemPtr problem, CSGraphPtr cs_graph, SSGraphMap ss_graphs, QPCCOptionsPtr options);
};

}