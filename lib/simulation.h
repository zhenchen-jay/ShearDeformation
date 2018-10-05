#ifndef SHEARDEFORMATION_SIMULATION_H
#define SHEARDEFORMATION_SIMULATION_H
#include <Eigen/Dense>
#include <vector>
#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"

class Simulation{
public:
    Simulation() : _is_initialized(false) {};

    ~Simulation() {};


public:
    void compute_deformed_surface(Eigen::MatrixXd &V0, Eigen::MatrixXi &F0, double nu, double E, double t);
    /*
     This function is aimed to compute the deformed surface based on the interior energy and external energy
     defined in this class.
     f = E_M + \lambda * E_external
     @param[in] V0: The vertice of the surface
     @param[in] F0: The triangle meshes
     @oaram[in] nu: The Passion ratio
     @param[in] E:  The material parameter
     @param[in] t:  Thickness

     @param[out] V0: The vertice of the surface after deformation
     @param[out] F0: The triangle meshes of the surface after deformation
     */

    void energy_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df);
    /*
     This function calculate the overall engergy and its gradient at the given position x.
     @param[in] x:     The variable
     @param[in] ptr:   A void pointer, which is required by alglib
     
     @param[out] f:    The overall energy function
     @param[out] df:   The gradient of the overall energy function
     
     */
    
private:

    // Membrane Energy

    double membrane_energy(Eigen::MatrixXd V);

    /*
     This function calculates the membrane energy
     @param[in] V: The vertice of the surface
     @return:      The Membrane Energy
     */
    
    double membrane_energy_default(Eigen::MatrixXd V);

    Eigen::VectorXd diff_membrane_energy(Eigen::MatrixXd V);

    /*
     This function calculates the direvative of the membrane energy
     @param[in] V0:     The vertice of the surface
     @param[out] dE:    The derivation of the Energy function
     */

    // Energy component cause by fixed points
    
     Eigen::VectorXd diff_membrane_energy_default(Eigen::MatrixXd V);

    double fixed_points_energy(Eigen::MatrixXd V);
    /*
     This function calculates the membrane energy
     @param[in] V: The vertice of the surface
     @return:      The Energy caused by the change of fixed point
     */

    Eigen::VectorXd diff_fixed_points_energy(Eigen::MatrixXd V);
    /*
     This function calculates the direvative of the membrane energy
     @param[in] V0:     The vertice of the surface
     @param[out] dE:    The derivation of the Energy function
     */
    
    // External Energy, like gravity
    
    double external_energy(Eigen::MatrixXd V);
    /*
     This function calculates the membrane energy
     @param[in] V: The vertice of the surface
     @return:      The energy caused by the external force
     */
    
    Eigen::VectorXd diff_external_energy(Eigen::MatrixXd V);
    /*
     This function calculates the direvative of the membrane energy
     @param[in] V0:     The vertice of the surface
     @param[out] dE:    The derivation of the Energy function
     */
    
public:
    void test_function(Eigen::VectorXd v, Eigen::MatrixXd V, double eps);
    void test(Eigen::MatrixXd &V, Eigen::MatrixXi F,double nu, double E, double t,double ratio);
    void initialize(Eigen::MatrixXd V0, Eigen::MatrixXi F0, double nu, double E, double t);
    

private:
    

    void calculate_fundamental_form1(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Matrix2d &I);

    /*
     This function calculates the corresponding matrix of the first fundamental form with the parametric space span{(0,0),(1,0) (0,1)}
     @param[in] {V0,V1,V2}: the vertice of the triangle
     @param[out] I: the corresponding matrix of the first fundamental form
     */

    double calculate_QI(Eigen::Matrix2d I_U, Eigen::Matrix2d I_D, Eigen::Vector2d v);

    /*
     This function calculates the quadratic form Q1 based on first fundamental form: Q1 = sqrt(I_U).inverse()*(I_D-I_U)sqrt(I_U).inverse(). Then it returns the Q(v).
     @param[in] I_U: The corresponding matrix for the first fundamental form of the undeformed surface.
     @param[in] I_D: The corresponding matrix for the first fundamental form of the deformed surface.
     @param[in] v:   The tangent vector.

     @return:        The value of QI(v)
     */

    Eigen::Vector3d diff_QI(Eigen::Matrix2d I_U, Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2, int flag,
                            Eigen::Vector2d e);

    /*
     This function calculates the derivative of the quadratic form Q1 based on first fundamental form: Q1 = sqrt(I_U).inverse()*(I_D-I_U)sqrt(I_U).inverse(), wrt. the vector v.
     @param[in] I_U: The corresponding matrix for the first fundamental form of the undeformed surface.
     @param[in] {v0,v1,v2}: The vertice of the triangle, which are used to determine dI_D
     @param[in] flag: The variable v_flag which belongs to {v0,v1,v2}
     @param[in] e:   The tangent vector.

     @return dQ: The value of dQ(v)
     */

    void calculate_coeff_mat(Eigen::Vector2d v0, Eigen::Vector2d v1, Eigen::Vector2d v2, Eigen::Matrix3d &q, double nu);

    /*
     This function calculates the coefficient matrix based on the base v0,v1,v2. This base when forming as a triangle should be anticlockeise
     */

    double
    get_coeff_mat_components(Eigen::Vector2d vj, Eigen::Vector2d vk, Eigen::Vector2d vp, Eigen::Vector2d vq, double nu);

    /*
     This function returns the value q(i,o), where q is the coefficient matrix
     @param[in] {vj,vk,vp,vq}: The dual edge of {ej,ek,ep,eq}, namely the edges rotated clockwise by 90 degrees.
     @param[in] nu: The Posion ratio.
     @return: q(i,o)
     */

    void
    diff_fundamental_form1(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2, int flag, Eigen::Matrix2d &dI_x,
                           Eigen::Matrix2d &dI_y, Eigen::Matrix2d &dI_z);
    /*
     This function calculates the derivative of the first fundamental form determined by v0,v1,v2 wrt. the vertex v_flag. In this function, we require that the flag belongs to {0,1,2}
     @param[in] {v0,v1,v2}: The vertice of the triangle.
     @param[in] flag: The variable v_flag which belongs to {v0,v1,v2}
     @param[out] {dI_x, dI_y, dI_z}: The derivative matrix wrt. x,y,z-coordinates.
     */

private:
    Eigen::MatrixXd V_U; // Vertice of undeformed surface
    Eigen::MatrixXd V_U_record; // Record the vertices infromation for the first initialization
    Eigen::MatrixXi F;  // Face information
    std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
    std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
    std::vector<Eigen::Matrix2d> IU_list; // list of the first fundamental form for the undeformed surface.
    std::vector<double> dA_list;            // list of the area ratio, which equals to sqrt(det(I_U))

    bool _is_initialized;

    // The base vectors
    Eigen::Vector2d e0;
    Eigen::Vector2d e1;
    Eigen::Vector2d e2;

    // Fixed points
    std::vector<Eigen::Vector3d> p_fixed;
    std::vector<int> p_fixed_index;

    Eigen::Matrix3d _q;  // The coefficient matrix q
    double _nu;  // The Possion Ratio
    double _A;  // Area of the parametric triangle
    double _alpha; // The coeffient of the Membrane energy
    double _beta;   // The coeffient of the Membrane energy
    double _t;  // The thickness
};


#endif //SHEARDEFORMATION_SIMULATION_H
