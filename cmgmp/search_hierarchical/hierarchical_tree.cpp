#include "hierarchical_tree.h"


HTree::HTree(const double& a_weight, const double& p_weight){
    this->angle_weight = a_weight;
    this->translation_weight = p_weight;
}

double HTree::dist(const Vector7d& q1, const Vector7d& q2){
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

int HTree::nearest_neighbor(const Vector7d& q)
{
    int near_idx;
    double min_d= DBL_MAX;
    for (int i = 0; i < nodes.size(); i++)
    {
    if (!nodes[i].is_alive){
            continue;
    }
    double d = this->dist(nodes[i].config, q);
    if (d < min_d)
    {
        near_idx = i;
        min_d = d;
    }
    }
    return near_idx;
}

int HTree::nearest_unextended_to_goal(const Vector7d& q){
    int near_idx;
    double min_d= DBL_MAX;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (!nodes[i].is_alive || nodes[i].is_extended_to_goal){
          continue;
        }
        double d = this->dist(nodes[i].config, q);
        if (d < min_d)
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

void HTree::backtrack(int last_node_idx, std::vector<int>* node_path)
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

void HTree::initial_node(HNode* n)
{
    n->parent = -1;
    nodes.push_back(*n);
    return;
}

void HTree::add_node(HNode* n, int parent_idx, HEdge* e)
{
    int node_idx = nodes.size();
    int edge_idx = edges.size();
    n->parent = parent_idx;
    n->edge = edge_idx;

    nodes.push_back(*n);
    edges.push_back(*e);
    return;
}


