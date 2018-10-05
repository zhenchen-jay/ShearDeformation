#ifndef MeshConnection_h
#define MeshConnection_h
#include <Eigen/Dense>
#include <vector>
namespace MeshConnection
{
    void vertex_triangle_adjacency(int vert_num, Eigen::MatrixXi F, std::vector<std::vector<int> > &VF, std::vector<std::vector<int> > &VFi);
    
    /*
    This function calculates the adjacent faces of the vertices using the function in libigl
    */
    
    
    void triangle_triangle_adjacency(Eigen::MatrixXi F, Eigen::MatrixXi &TT, Eigen::MatrixXi &TTi);
    
    /*
    This function calculates the adjacent faces of the face using the function in libigl 
    */

   void edge_triangle_adjacency(Eigen::MatrixXd V, Eigen::MatrixXi F);
   /*
   This function calculates the adjacent faces of the edges in a given mesh(V,F)
   */

    
    
}
#endif
