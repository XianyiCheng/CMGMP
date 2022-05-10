#include "contact_mode_enumeration.h"

#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
    #include "contact_constraints.h"
#endif

#include <modus/modes/enumerate.hpp>
#include <modus/modes/geometry/interior_point.hpp>
#include <modus/common/linear_algebra.hpp>
// #include <modus/common/serialization.hpp>
// #include <modus/common/logging.hpp>
#include <modus/modes/geometry/incidence_graph.hpp>

void cs_mode_enumeration(const MatrixXd& A, std::vector<VectorXi>* cs_modes){
   
    if (A.rows() == 1){
        VectorXi m1(1);
        m1[0] = 1;
        VectorXi m0(1);
        m0[0] = 0;
        cs_modes->push_back(m0);
        cs_modes->push_back(m1);
        return;
    }
    ModeEnumerationOptions* options = new ModeEnumerationOptions();
    IncidenceGraph* graph = 0;
    double eps = 1e-4;
    graph = modus::EnumerateCSModes(-A, eps, options);
    if (graph == 0){
        return;
    }
    // enumerate_cs_modes(A, b, 1e-8, options);
    std::vector<std::string> sv = graph->get_proper_sign_vectors();

    for(auto& s: sv){
        VectorXi ss(s.size());
        for(int i = 0; i < s.size(); i++){
            if(s[i] == '0'){
                ss[i] = 0;
            } else {
                ss[i] = 1;
            }
        }
        cs_modes->push_back(ss);
    }
    delete graph;
}

void ss_mode_enumeration(const MatrixXd& A, const MatrixXd& T, VectorXi cs_mode, std::vector<VectorXi>* ss_modes){
    
    ModeEnumerationOptions* options = new ModeEnumerationOptions();

    double eps = 1e-4;

    std::string s;
    for(int i = 0; i < cs_mode.size(); i++){
        if(cs_mode[i] == 0){
            s.push_back('0');
        } else {
            s.push_back('-');
        }
    }

    

    if (cs_mode.sum() == cs_mode.size()){
        VectorXi ss(s.size()*3);
        ss.setZero();
        ss.block(0,0,s.size(),1) = cs_mode;
        ss_modes->push_back(ss);
        return;
    }

    IncidenceGraph* ss_graph = 0;
    // std::cout << "A" << A << "\n T" << T << std::endl;
    ss_graph = modus::EnumerateSSModes(-A, -T, s, eps, options);
    if (ss_graph == 0){
        std::cout << "FAIL TO ENUMERATE SS MODES! CS MODE: " << cs_mode.transpose() << std::endl;
        return;
    }
    // std::cout <<"ss mode:" << std::endl;

    for (int m = 0; m < ss_graph->rank(); m++) {
        for (int n = 0; n < ss_graph->rank(m).size(); n++) {
            Node* uu = ss_graph->rank(m)[n];
            // std::cout << uu->sign_vector << std::endl;
            std::string uu_s = uu->sign_vector;
            std::string s_uu = uu_s.substr(0,s.size());
            if (s.compare(s_uu)!=0){
                continue;
            }
            VectorXi ss(s.size()*3);
            ss.setZero();
            ss.block(0,0,s.size(),1) = cs_mode;
            int counter = 0;
            for(int ss_i = s.size(); ss_i < 3*s.size(); ss_i++){
                if(s[(ss_i-s.size())/2] != '0'){
                    continue;
                } else if(uu_s[s.size()+counter] == '+') {
                    ss[ss_i] = -1;
                } else if(uu_s[s.size()+counter] == '-') {
                    ss[ss_i] = 1;
                } else {
                    ss[ss_i] = 0;
                }
                counter++;
            }
            ss_modes->push_back(ss);
            // std::cout << s << " " << uu_s << " " << ss.transpose() << std::endl;
        }
        
    }
    delete ss_graph;
    return;
}

void all_mode_enumeration(const MatrixXd& A, const MatrixXd& T, std::vector<VectorXi>* all_modes){
    // std::cout << -A << std::endl;
    // std::cout << -T << std::endl;

    
    
    if (A.rows() == 1){

        VectorXi m0(1);
        m0[0] = 0;

        ss_mode_enumeration(A, T, m0, all_modes);

        VectorXi m1(3);
        m1 << 1,0,0;

        all_modes->push_back(m1);

        return;

    }

    
    ModeEnumerationOptions* options = new ModeEnumerationOptions();
    IncidenceGraph* graph = 0;
    double eps = 1e-4;
    graph = modus::EnumerateCSModes(-A, eps, options);
    if (graph == 0){
        return;
    }
    // enumerate_cs_modes(A, b, 1e-8, options);
    std::vector<std::string> sv = graph->get_proper_sign_vectors();

    // Check interior points against sign vectors.
    for (int k = 0; k < graph->rank(); k++) {
        // std::cout << "rank " << k << std::endl;
        for (int j = 0; j < graph->rank(k).size(); j++) {
            Node* u = graph->rank(k)[j];
            // std::string sv = get_sign_vector(A * u->interior_point, 1e-8);
            // std::cout << sv << std::endl;
            // std::cout <<"cs mode:"<< u->sign_vector << std::endl;
            // std::cout << *u << std::endl << std::endl;
            // std::cout << A << std::endl;
            // std::cout << A * u->interior_point - b << std::endl << std::endl;

            std::string s = u->sign_vector;
            VectorXi cs(s.size());
            for(int cs_i = 0; cs_i < s.size(); cs_i++){
                if(s[cs_i] == '0'){
                    cs[cs_i] = 0;
                } else {
                    cs[cs_i] = 1;
                }
            }

            if (cs.sum() == s.size()){
                VectorXi ss(s.size()*3);
                ss.setZero();
                ss.block(0,0,s.size(),1) = cs;
                all_modes->push_back(ss);
                continue;
            }

            IncidenceGraph* ss_graph = 0;

            
            ss_graph = modus::EnumerateSSModes(-A, -T, u->sign_vector, eps, options);
            // std::cout <<"ss mode:" << std::endl;

            if (ss_graph == 0){
                continue;
            }
        
            for (int m = 0; m < ss_graph->rank(); m++) {
                for (int n = 0; n < ss_graph->rank(m).size(); n++) {
                    Node* uu = ss_graph->rank(m)[n];
                    // std::cout << uu->sign_vector << std::endl;
                    std::string uu_s = uu->sign_vector;
                    std::string s_uu = uu_s.substr(0,s.size());
                    if (s.compare(s_uu)!=0){
                        continue;
                    }
                    VectorXi ss(s.size()*3);
                    ss.setZero();
                    ss.block(0,0,s.size(),1) = cs;
                    int counter = 0;
                    for(int ss_i = s.size(); ss_i < 3*s.size(); ss_i++){
                        if(s[(ss_i-s.size())/2] != '0'){
                            continue;
                        } else if(uu_s[s.size()+counter] == '+') {
                            ss[ss_i] = -1;
                        } else if(uu_s[s.size()+counter] == '-') {
                            ss[ss_i] = 1;
                        } else {
                            ss[ss_i] = 0;
                        }
                        counter++;
                    }
                    all_modes->push_back(ss);
                    // std::cout << s << " " << uu_s << " " << ss.transpose() << std::endl;
                }
                
            }

            delete ss_graph;
        }
    }
    delete graph;
}

