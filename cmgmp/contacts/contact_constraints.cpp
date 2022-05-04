#include "contact_constraints.h"

#ifndef CONTACTKINEMATICS_H
#define CONTACTKINEMATICS_H
    #include "contact_kinematics.h"
#endif

#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#ifndef SAMPLE_H
#define SAMPLE_H
    #include "../utilities/sample.h"
#endif

#define PI 3.1415926

void numeric_stablize(MatrixXd& M){
    for(int i=0; i<M.rows(); i++){
        for(int j = 0; j<M.cols(); j++){
            if(abs(M(i,j)) < 1e-10){
                M(i,j)=0;
            }
        }
    }
}

FrictionConeLinearization::FrictionConeLinearization(int n) 
{
    double thr = 1e-3;
    number_of_sliding_planes = n;
    D.resize(n, Eigen::NoChange);
    D.setZero();
    D_dual.resize(2*n, Eigen::NoChange);
    D_dual.setZero();
    cones.resize(2*2*n, Eigen::NoChange);
    vectors.resize(2*n, Eigen::NoChange);
    cone_modes.resize(2*n, n);
    vector_modes.resize(2*n, n);

    for (int i = 0; i < n; i++){
        D(i,0) = cos(PI*i/n);
        D(i,1) = sin(PI*i/n);
    }

    for (int j = 0; j < 2*n; j++){

        D_dual(j,0) = cos(PI*j/n + PI/2);
        D_dual(j,1) = sin(PI*j/n + PI/2);

        vectors.row(j) = D_dual.row(j);

        MatrixXd cone(2,6);

        if (j == 0) {
            cone.setZero();
            cone(0,0) = cos(PI*(2*n-1)/n + PI/2);
            cone(0,1) = sin(PI*(2*n-1)/n + PI/2);
            cone(1,0) = cos(PI*j/n + PI/2);
            cone(1,1) = sin(PI*j/n + PI/2);
        } else {
            cone = D_dual.block<2,6>(j-1,0);
        }
        cones.block<2,6>(2*j,0) = cone;

        MatrixXd cone_mode = (cone.colwise().sum()/2)*D.transpose();
        MatrixXd vector_mode = D_dual.row(j)*D.transpose();

        for (int k = 0; k < n; k++){

            if (cone_mode(0,k) > thr){
                cone_modes(j,k) = 1;
            } else if (cone_mode(0,k) < -thr){
                cone_modes(j,k) = -1;
            } else {
                cone_modes(j,k) = 0;
            }
            
            if (vector_mode(0,k) > thr){
                vector_modes(j,k) = 1;
            } else if (vector_mode(0,k) < -thr){
                vector_modes(j,k) = -1;
            } else {
                vector_modes(j,k) = 0;
            }
        }
    }

}

MatrixXd FrictionConeLinearization::getVectorsbyMode(const VectorXi& mode){


    bool is_vector = false;
    bool is_cone = false;
    int index = -1;
    MatrixXd vs;
    for (int i = 0; i < 2*number_of_sliding_planes; i++){
        is_vector = true;
        is_cone = true;
        for (int j = 0; j < number_of_sliding_planes; j++){
            if (vector_modes(i,j) != mode(j)){
                is_vector = false;
            }
            if (cone_modes(i,j) != mode(j)){
                is_cone = false;
            }
        }
        if (is_vector || is_cone){
            vs.resize(2,6);
            vs.setZero();
            index = i;
            break;
        }
    }
    if (is_cone){
        vs = cones.block(2*index,0,2,6);
    }
    if (is_vector){
        vs.row(0) = vectors.row(index);
    }
    return vs;
}

ContactConstraints::ContactConstraints (int num_sliding_planes) {
    
    friction_cone = new FrictionConeLinearization(num_sliding_planes);
    
    basis.setZero();
    for (int i = 0; i < 3; i++){
        basis(i,i) = 1;
    }
}

void ContactConstraints::ModeConstraints(const std::vector<ContactPoint>&pts ,VectorXi mode, double mu, Vector6d f_external, 
                                    MatrixXd* A_pr, VectorXd* b_pr, MatrixXd* G_pr, VectorXd* h_pr)
{
    int n_pts = pts.size();
    const int n = this->friction_cone->number_of_sliding_planes;

    int n_var = 6;
    for (int i = 0; i < n_pts; i++){
        int cs_mode = mode(i);
        if (cs_mode == 0){
            n_var += 3;
        }
    }

    MatrixXd A_static(6, n_var - 6);

    MatrixXd A(3*n_pts+6, n_var);
    MatrixXd G((2+n*2+1)*n_pts, n_var);
    A_static.setZero();
    A.setZero();
    G.setZero();

    int counter_A = 0;
    int counter_G = 0;
    int k = 0; // counter for contacting contacts

    for (int i = 0; i < n_pts; i++){

        int cs_mode = mode(i);
        MatrixXi ss_mode(n,1);
        ss_mode = mode.block(n_pts + i*n, 0,n,1);

        Matrix6d Adgco = contact_jacobian(pts[i].p, pts[i].n);

        // std::cout << "Adgco\n" << Adgco << std::endl;
        
        if (cs_mode == 1){ // separate

            G.block(counter_G,0,1,6) = basis.row(2)*Adgco;
            counter_G += 1;

        } else if (ss_mode.cwiseAbs().sum() == 0) // all fixed
        {
            A.block(counter_A, 0, 3, 6) = basis*Adgco;
            counter_A += 3;

            //TODO: use pyramind approximation for now, fix with the actual approximation later
            MatrixXd pyramid(4,3);
            pyramid << 1,1,mu,
                       -1,1,mu,
                       1,-1,mu,
                       -1,-1,mu;
                    
            G(counter_G, 6+k*3 + 2) = 1;
            G.block(counter_G+1, 6+k*3, 4, 3) = pyramid;
            counter_G += 1+2*2;

            A_static.block(0, 3*k, 6, 3) = (basis*Adgco).transpose();
            ++k;

        } else { // sliding

            // normal velocity constraints n_v = 0
            A.block(counter_A, 0, 1, 6) = basis.row(2) * Adgco;
            counter_A += 1;

            // sliding velocity constraints
            MatrixXd Sv(n, 6);
            Sv = friction_cone->D * Adgco;

            for (int j = 0; j < n; j++){
                int m = ss_mode(j);
                if (m == 0){
                    A.block(counter_A, 0, 1, 6) = Sv.row(j);
                    counter_A += 1;
                } else if (m == 1) {
                    G.block(counter_G, 0, 1, 6) = Sv.row(j);
                    counter_G += 1;
                } else { // m == -1
                    G.block(counter_G, 0, 1, 6) = -Sv.row(j);
                    counter_G += 1;
                }
            }

            // friction force
            //
            /*
            Eigen::Matrix<double, 5, 3> fc;
            fc.setZero();
            fc.block(0,0,3,3) = Matrix3d::Identity();
            fc.block(3,0,2,2) = -Eigen::Matrix2d::Identity();
            fc(3,2) = mu;
            if (ss_mode.cwiseAbs().minCoeff() > 0){
                fc(4,2) = mu;
            }
            */

            // fix sliding force bug
            Eigen::Matrix<double, 5, 3> fc;
            
            if (ss_mode.cwiseAbs().minCoeff() > 0){
                fc << 1,0,0,
                      0,1,0,
                      0,0,1,
                      -1,-1,1.01*mu,
                      1,1,-0.99*mu;
            } else {
                fc << 1,0,0,
                      0,1,0,
                      0,0,1,
                      -1,0,1.01*mu,
                      1,0,-0.99*mu;
            }

            G.block(counter_G, 6+k*3, 5, 3) = fc;
            counter_G += 5;

            MatrixXd modecone = friction_cone->getVectorsbyMode(ss_mode);
            Eigen::Matrix<double, 3, 6> force_basis;
            force_basis.block(0,0,2,6) = -modecone;
            force_basis.block(2,0,1,6) = basis.row(2);

            A_static.block(0, k*3, 6, 3) = (force_basis*Adgco).transpose();

            ++k;
        }
    }

    A.block(counter_A, 6, 6, n_var - 6) = A_static;
    counter_A += 6;

    A.conservativeResize(counter_A, n_var);
    G.conservativeResize(counter_G, n_var);
    
    VectorXd b(counter_A);
    VectorXd h(counter_G);
    b.setZero();
    b.block(counter_A-6, 0, 6, 1) = -f_external;
    h.setZero();

    *A_pr = A;
    *b_pr = b;
    *G_pr = G;
    *h_pr = h;

    return;
}

void deleteZeroRows(const MatrixXd& A, const VectorXd& b, MatrixXd* A_pr, VectorXd* b_pr){
    
    MatrixXd AA(A.rows(), A.cols());
    VectorXd bb(b.rows());

    AA.setZero();
    bb.setZero();

    int k = 0;

    bool isnonzero = false;
    for (int i = 0; i < A.rows(); i++){

        isnonzero = false;
        for (int j = 0; j < A.cols(); j++){
            if (A(i,j)!=0){
                isnonzero = true;
                break;
            }
        }

        if (isnonzero){
            AA.row(k) = A.row(i);
            bb(k) = b(i);
            k += 1;
        }
    }

    AA.conservativeResize(k, Eigen::NoChange);
    bb.conservativeResize(k);

    *A_pr = AA;
    *b_pr = bb;

}

void mergeDependentRows(const MatrixXd& A, const VectorXd& b, MatrixXd* A_pr, VectorXd* b_pr){
    
    MatrixXd AA(A.rows(), A.cols());
    VectorXd bb(b.rows());

    AA.setZero();
    bb.setZero();

    int k = 0;

    bool isdependent = false;
    for (int i = 0; i < A.rows(); i++){

        isdependent = false;

        MatrixXd Ai(1,A.cols());
        Ai = A.block(i,0,1,A.cols());
        Ai = Ai/Ai.norm();
        MatrixXd M(k+1, A.cols());
        M = AA.block(0,0,k+1, A.cols());
        for (int ii = 0; ii < k+1; ii++){
            double sumabs = (Ai - M.row(ii)/M.row(ii).norm()).cwiseAbs().sum();
            if (sumabs < 1e-5){
                isdependent = true;
                break;
            }
        }


        if (!isdependent){
            AA.row(k) = A.row(i);
            bb(k) = b(i);
            k += 1;
        }
    }

    AA.conservativeResize(k, Eigen::NoChange);
    bb.conservativeResize(k);

    *A_pr = AA;
    *b_pr = bb;

}

void mergeManipulatorandEnvironmentConstraints(const MatrixXd& A_mnp, const VectorXd& b_mnp, const MatrixXd& G_mnp, const VectorXd& h_mnp,
    const MatrixXd& A_env, const VectorXd& b_env, const MatrixXd& G_env, const VectorXd& h_env,
    MatrixXd* A_pr, VectorXd* b_pr, MatrixXd* G_pr, VectorXd* h_pr)
{


    int n_var = A_mnp.cols() + A_env.cols() - 6;
    int A_rows = A_mnp.rows() + A_env.rows() - 6;
    int G_rows = G_mnp.rows() + G_env.rows();

    MatrixXd A(A_rows, n_var);
    VectorXd b(A_rows);
    MatrixXd G(G_rows, n_var);
    VectorXd h(G_rows);
        
    A.setZero();
    G.setZero();
    b.setZero();
    h.setZero();


    // velocity part
    A.block(0,0,A_env.rows() - 6, A_env.cols()) = A_env.block(0, 0, A_env.rows() - 6, A_env.cols());
    A.block(A_env.rows() - 6, A_env.cols() , A_mnp.rows() - 6, A_mnp.cols() - 6) = A_mnp.block(0, 6, A_mnp.rows() - 6, A_mnp.cols() - 6);
    
    // force balance part
    A.block(A_rows - 6, 0, 6, A_env.cols()) = A_env.block(A_env.rows()-6, 0, 6, A_env.cols());
    A.block(A_rows - 6, A_env.cols(), 6, A_mnp.cols()-6) = A_mnp.block(A_mnp.rows()-6, 6, 6, A_mnp.cols()-6);

    b.block(A_rows - 6, 0, 6, 1) = b_env.block(b_env.rows() - 6, 0, 6, 1);

    G.block(0,0,G_env.rows(), G_env.cols()) = G_env;
    G.block(G_env.rows(), G_env.cols() , G_mnp.rows(), G_mnp.cols() - 6) = G_mnp.block(0, 6, G_mnp.rows(), G_mnp.cols() - 6);

    deleteZeroRows(A, b, A_pr, b_pr);
    deleteZeroRows(G, h, G_pr, h_pr);

    mergeDependentRows(*A_pr, *b_pr, A_pr, b_pr);
   
    // *A_pr = A;
    // *b_pr = b;
    // *G_pr = G;
    // *h_pr = h;
    
    return;
}

void mergeManipulatorandEnvironmentConstraints_relax(const MatrixXd& A_mnp, const VectorXd& b_mnp, const MatrixXd& G_mnp, const VectorXd& h_mnp,
    const MatrixXd& A_env, const VectorXd& b_env, const MatrixXd& G_env, const VectorXd& h_env,
    MatrixXd* A_pr, VectorXd* b_pr, MatrixXd* G_pr, VectorXd* h_pr)
{


    int n_var = A_mnp.cols() + A_env.cols() - 6;
    int A_rows = A_mnp.rows() + A_env.rows() - 6 - 6;
    int G_rows = G_mnp.rows() + G_env.rows() + 12;

    MatrixXd A(A_rows, n_var);
    VectorXd b(A_rows);
    MatrixXd G(G_rows, n_var);
    VectorXd h(G_rows);

    MatrixXd A_s(6, n_var);
    VectorXd b_s(6);
        
    A.setZero();
    G.setZero();
    b.setZero();
    h.setZero();


    // velocity part
    A.block(0,0,A_env.rows() - 6, A_env.cols()) = A_env.block(0, 0, A_env.rows() - 6, A_env.cols());
    A.block(A_env.rows() - 6, A_env.cols() , A_mnp.rows() - 6, A_mnp.cols() - 6) = A_mnp.block(0, 6, A_mnp.rows() - 6, A_mnp.cols() - 6);
    
    // force balance part
    A_s.block(0, 0, 6, A_env.cols()) = A_env.block(A_env.rows()-6, 0, 6, A_env.cols());
    A_s.block(0, A_env.cols(), 6, A_mnp.cols()-6) = A_mnp.block(A_mnp.rows()-6, 6, 6, A_mnp.cols()-6);

    b_s = b_env.block(b_env.rows() - 6, 0, 6, 1);

    G.block(0,0,G_env.rows(), G_env.cols()) = G_env;
    G.block(G_env.rows(), G_env.cols() , G_mnp.rows(), G_mnp.cols() - 6) = G_mnp.block(0, 6, G_mnp.rows(), G_mnp.cols() - 6);

    VectorXd sigma = VectorXd::Constant(b_s.size(), 1e-4);
    G.block(G_env.rows()+G_mnp.rows(), 0 , 6, n_var) = A_s;
    G.block(G_env.rows()+G_mnp.rows()+6, 0 , 6, n_var) = -A_s;
    h.block(G_env.rows()+G_mnp.rows(), 0 , 6, 1) = b_s - sigma;
    h.block(G_env.rows()+G_mnp.rows()+6, 0 , 6,1) = - (b_s + sigma);

    deleteZeroRows(A, b, A_pr, b_pr);
    deleteZeroRows(G, h, G_pr, h_pr);

    mergeDependentRows(*A_pr, *b_pr, A_pr, b_pr);
   
    // *A_pr = A;
    // *b_pr = b;
    // *G_pr = G;
    // *h_pr = h;
    
    return;
}

void ContactConstraints::NormalVelocityConstraints(const std::vector<ContactPoint>& pts, MatrixXd* A_pr, VectorXd* b_pr){
    //  Create halfspace inequalities, Aq̇ - b ≤ 0

    int n_pts = pts.size();
    MatrixXd A(n_pts, 6);
    VectorXd b(n_pts);
    A.setZero();
    b.setZero();
    
    for(int i = 0; i < n_pts; i++){
        Matrix6d Adgco = contact_jacobian(pts[i].p, pts[i].n);
        A.row(i) = this->basis.row(2)*Adgco;
    }
    
    A = -A;

    // numeric_stablize(A);

    *A_pr = A;
    *b_pr = b;
}

void ContactConstraints::TangentVelocityConstraints(const std::vector<ContactPoint>& pts, MatrixXd* T_pr, VectorXd* t_pr){
    
    int n_pts = pts.size();
    int n_slide = this->friction_cone->number_of_sliding_planes;
    MatrixXd T(n_pts*n_slide, 6);
    VectorXd t(n_pts*n_slide);
    T.setZero();
    t.setZero();
    
    for(int i = 0; i < n_pts; i++){
        Matrix6d Adgco = contact_jacobian(pts[i].p, pts[i].n);
        T.block(i*n_slide,0,n_slide,6) = this->friction_cone->D*Adgco;
    }
    
    T = -T;

    // numeric_stablize(T);

    *T_pr = T;
    *t_pr = t;
}

void copy_points(const std::vector<ContactPoint>& pts, std::vector<ContactPoint>* pts_new){
    for(auto& pt:pts){
        pts_new->push_back(pt);
    }
}

void print_contacts(const std::vector<ContactPoint>& pts){
    std::cout << "Contact points: " << std::endl;
    for (auto cp: pts){
        std::cout << "point:" << cp.p.transpose() 
                << ", normal:" << cp.n.transpose() 
                << ", depth: " << cp.d << std::endl;
    }
}
