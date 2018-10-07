//
//  GeometryFeature.h
//  ShearDeformation_bin
//
//  Created by Zhen Chen on 10/5/18.
//

#ifndef GeometryFeature_h
#define GeometryFeature_h
#include<Eigen/Dense>
#include<vector>
namespace GeoFeature
{
    void face_normal(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int start_ind, Eigen::Vector3d &nf, std::vector<Eigen::Vector3d> &dn);
    
    /*
     This function calculates the face normals and its derivative w.r.t the given vertex, which can be decided by the start_ind. That is, start_ind = i, then V = Vi
     @param[in] {V0, V1, V2}: The vertices of the triangle.
     @param[in] start_ind: The flag of variable V which belongs to {V0,V1,V2}, hat is, start_ind = i, then V = Vi. If start_ind = -1, calculate the normal and return dn as zero
     
     @param[out] nf: The face normal.
     @param[out] dn: The derivative of this normal, w.r.t. the vertex decided by start_ind
     
     */
    
    
    void calculate_first_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Matrix2d &I);
    
    /*
     This function calculates the corresponding matrix of the first fundamental form with the parametric space span{(0,0),(1,0) (0,1)}
     @param[in] {V0,V1,V2}: the vertices of the triangle
     @param[out] I: the corresponding matrix of the first fundamental form
     */
    
    void
    diff_first_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int start_ind, std::vector<Eigen::Matrix2d> &dI);
    /*
     This function calculates the derivative of the first fundamental form determined by V0,V1,V2 w.r.t. the vertex v_flag. In this function, we require that the flag belongs to {0,1,2}
     @param[in] {V0,V1,V2}: The vertices of the triangle.
     @param[in] start_ind: The flag of variable V which belongs to {V0,V1,V2}, that is, start_ind = i, then V = Vi
     @param[out] dI: The vector fo the derivative matrix wrt. x,y,z-coordinates.
     */
    
    void calculate_second_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Vector3d V_0, Eigen::Vector3d V_1, Eigen::Vector3d V_2, std::vector<bool> real_pts, Eigen::Matrix2d &II);
    /*
     This function calculates the corresponding matrix of the second fundamental form with the parametric space span{(0,0),(1,0) (0,1)}
    @param[in] {V_0,V_1,V_2}: The vertices of the adjacent triangles, s.t. {V_0,V2,V1},{V_1,V0,V2},{V_2,V1,V0} are adjacent triangles on the given triangle mesh.
     @param[in] {n0,n1,n2}: the corresponding edge normal of the triangle
     @param[in] real_pts: check whether V_k is existed, which will happen on the boundary
     
     @param[out] II: the corresponding matrix of the first fundamental form
     */
    
    void diff_second_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Vector3d V_0, Eigen::Vector3d V_1, Eigen::Vector3d V_2, int start_ind, std::vector<bool> real_pts, std::vector<Eigen::Matrix2d> &dII);
    /*
     This function calcualtes the gradient of the second fundamental form w.r.t. the vertex which lies on this triangle.
     @param[in] {V0,V1,V2}: The vertices of the triangle
     @param[in] {V_0,V_1,V_2}: The vertices of the adjacent triangles, s.t. {V_0,V2,V1},{V_1,V0,V2},{V_2,V1,V0} are adjacent triangles on the given triangle mesh.
     @param[in] start_ind: The flag of variable V which belongs to {V0,V1,V2} or {V_0,V_1,V_2}, that is, start_ind = i(0<=i<=2), then V = Vi, otherwise, start_ind = i+3 (0<=i<=2), then V = V_i
     @param[in] real_pts: check whether V_k is existed, which will happen on the boundary
     
     @param[out] dII: The gradient of the second fundamental form w.r.t to the vertex V_{start_ind}
     
     */
    
    
    // Test functions which can be used to check whether the derivative is correct
    
    void test_face_normal();
    
    void test_second_fundamental_form();
    
}
#endif /* GeometryFeature_h */
