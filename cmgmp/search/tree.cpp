#include "tree.h"


Tree::Tree(const double& a_weight, const double& p_weight){
    this->angle_weight = a_weight;
    this->translation_weight = p_weight;
}

double Tree::dist(const Vector7d& q1, const Vector7d& q2){
    Vector3d p1(q1[0],q1[1],q1[2]);
    Vector3d p2(q2[0],q2[1],q2[2]);
    Quaterniond quat1(q1[6], q1[3], q1[4], q1[5]); // Quaterniond: w, x, y, z
    Quaterniond quat2(q2[6], q2[3], q2[4], q2[5]);
    Vector3d dp = p1 - p2;
    double dtrans = dp.norm();
    double angle = angBTquat(quat1, quat2); 
    
    double d = this->translation_weight*dtrans + this->angle_weight*angle;
    return d;

}

int Tree::nearest_neighbor(const Vector7d& q)
{
  int near_idx;
  double min_d= DBL_MAX;
  for (int i = 0; i < nodes.size(); i++)
  {
      double d = this->dist(nodes[i].config, q);
      if (d < min_d)
      {
          near_idx = i;
          min_d = d;
      }
  }
  return near_idx;
}

int Tree::nearest_unextended_to_goal(const Vector7d& q){
    int near_idx;
    double min_d= DBL_MAX;
    for (int i = 0; i < nodes.size(); i++)
    {
        double d = this->dist(nodes[i].config, q);
        if ((d < min_d) && (!nodes[i].is_extended_to_goal))
        {
            near_idx = i;
            min_d = d;
        }
    }
    if (min_d == DBL_MAX){
        return -1;
    } 
    
    return near_idx;
}

// void Tree::nearest_neighbors(Vector7d q[7], std::vector<int>* neighbors)
// {

//   double min_d= DBL_MAX;
//   for (int i = 0; i < nodes.size(); i++)
//   {
//       double d = dist(nodes[i].config, q);
//       if (d < min_d)
//       {
//           min_d = d;
//       }
//   }

//   for (int i = 0; i < nodes.size(); i++)
//   {

//       double d = dist(nodes[i].config, q);
//       if ((d - min_d) < 1e-5)
//       {
//           neighbors->push_back(i);
//       }
    
//     }
//   return;
// }

void Tree::backtrack(int last_node_idx, std::vector<int>* node_path)
{
    int cur_node = last_node_idx;
    node_path->push_back(cur_node);
    bool is_start_node = nodes[cur_node].parent != -1;
    while(is_start_node){
        cur_node = nodes[cur_node].parent;
        node_path->push_back(cur_node);
        is_start_node = nodes[cur_node].parent != -1;
    }
    return;
}

// void Tree::neighborhood(int node_idx, double radius, std::vector<int>* neighbors)
// {

//   for (int i = 0; i < nodes.size(); i++)
//   {
//       if (i==node_idx){continue;}
//       double d = dist(nodes[i].config, nodes[node_idx].config);
//       if (d <= radius)
//       {
//           neighbors->push_back(i);
//       }
    
//     }
//     return;
// }

void Tree::initial_node(Node* n)
{
    n->parent = -1;
    nodes.push_back(*n);
    return;
}

void Tree::add_node(Node* n, int parent_idx, Edge* e)
{
    int node_idx = nodes.size();
    int edge_idx = edges.size();
    n->parent = parent_idx;
    n->edge = edge_idx;

    // if (parent_idx!=-1){
    //     nodes[parent_idx].children.push_back(node_idx);
    // }

    nodes.push_back(*n);
    edges.push_back(*e);
    return;
}

// void Tree::remove_parent(int node_idx) {
//     int parent_idx = nodes[node_idx].parent;
//     for (std::vector<int>::iterator it = nodes[parent_idx].children.begin();
//             it != nodes[parent_idx].children.end(); ++it) {
//         if (node_idx == *it) {
//             nodes[parent_idx].children.erase(it);
//             break;
//         }
//     }
//     nodes[node_idx].parent = -1;
//     return;
// }

// void Tree::set_parent(int parent, int child) {
//     nodes[parent].children.push_back(child);
//     nodes[child].parent = parent;
// }
