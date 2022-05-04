
#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

void cs_mode_enumeration(const MatrixXd& A, std::vector<VectorXi>* cs_modes);

void all_mode_enumeration(const MatrixXd& A, const MatrixXd& T, std::vector<VectorXi>* all_modes);

void ss_mode_enumeration(const MatrixXd& A, const MatrixXd& T, VectorXi cs_mode, std::vector<VectorXi>* ss_modes);