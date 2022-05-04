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


struct Node {
    Vector7d config; // x, y, z, qx, qy, qz, qw
    int parent = -1; 
    int edge = -1; // index of edge from parent to this
    bool is_extended_to_goal = false;
    // std::vector<ContactPoint> fingertips;
    VectorXd manipulator_config;
    std::vector<ContactPoint> envs;
    std::vector<VectorXi> modes;
    // std::vector<int> children; // index of children nodes
    // double cost = 0;
    Node(Vector7d data, VectorXd mnp_config){
        config = data;
        manipulator_config = mnp_config;
      }
};

struct Edge {
  VectorXi mode;
  bool manipulator_collide;
  std::vector<Vector7d> path;
  Edge(VectorXi m, bool m_collide, std::vector<Vector7d>& p){
    mode = m;
    manipulator_collide = m_collide;
    path = p;
  }
};

class Tree{
public:
  std::vector<Node> nodes;
  std::vector<Edge> edges;
  double angle_weight;
  double translation_weight;
  Tree(const double& a_weight, const double& p_weight);
  double dist(const Vector7d& q1, const Vector7d& q2);
  int nearest_neighbor(const Vector7d& q);
  int nearest_unextended_to_goal(const Vector7d& q);
  // void nearest_neighbors(Vector7d q, std::vector<int>* nearest_neighbors);
  void backtrack(int last_node_idx, std::vector<int>* node_path);
  // void neighborhood(int node_idx, double radius, std::vector<int>* neighbors);
  void add_node(Node* n, int parent, Edge* e);
  void initial_node(Node* n);
  // void remove_parent(int node_idx);
  // void set_parent(int parent, int child);
};
