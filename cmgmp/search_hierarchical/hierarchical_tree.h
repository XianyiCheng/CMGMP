#pragma once
#include <vector>
#include <cfloat>
#include <cstdlib>

#ifndef UTILS_H
#define UTILS_H
  #include "../utilities/utilities.h"
#endif 

#ifndef CONTACTCONSTRAINTS_H
#define CONTACTCONSTRAINTS_H
  #include "../contacts/contact_constraints.h"
#endif

#define QUASISTATIC 0
#define QUASIDYNAMIC 1
#define DYNAMIC 2


class HNode {
    Vector7d config; // x, y, z, qx, qy, qz, qw
    int parent = -1; 
    int edge = -1; // index of edge from parent to this

    bool is_extended_to_goal = false;

    bool is_alive = true;

    int dynamic_type = QUASISTATIC;

    std::vector<ContactPoint> envs;
    std::vector<VectorXi> modes;

    HNode(Vector7d data){
        config = data;
      }
    
    HNode(Vector7d data, int type){
        config = data;
        dynamic_type = type;
      }
};

class HEdge {
  VectorXi mode;
  bool manipulator_collide;
  std::vector<Vector7d> path;
  HEdge(VectorXi m, bool m_collide, std::vector<Vector7d>& p){
    mode = m;
    manipulator_collide = m_collide;
    path = p;
  }
};

class HTree{
public:
  std::vector<HNode> nodes;
  std::vector<HEdge> edges;
  double angle_weight;
  double translation_weight;
  HTree(const double& a_weight, const double& p_weight);
  double dist(const Vector7d& q1, const Vector7d& q2);
  int nearest_neighbor(const Vector7d& q);
  int nearest_unextended_to_goal(const Vector7d& q);

  void backtrack(int last_node_idx, std::vector<int>* node_path);

  void add_node(Node* n, int parent, Edge* e);
  void initial_node(Node* n);

};
