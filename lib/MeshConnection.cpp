#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include "MeshConnection.h"

void MeshConnection::vertex_triangle_adjacency(int vert_num, Eigen::MatrixXi F, std::vector<std::vector<int> > &VF, std::vector<std::vector<int> > &VFi)
{
    igl::vertex_triangle_adjacency(vert_num,F,VF,VFi);
}

void MeshConnection::triangle_triangle_adjacency(Eigen::MatrixXi F, Eigen::MatrixXi &TT, Eigen::MatrixXi &TTi)
{
    igl::triangle_triangle_adjacency(F,TT,TTi);
}

void MeshConnection::edge_triangle_adjacency(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    
}

