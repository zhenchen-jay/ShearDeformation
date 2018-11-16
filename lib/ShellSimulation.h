//
// Created by Zhen Chen on 10/3/18.
//

#ifndef SHEARDEFORMATION_SHELLSIMULATION_H
#define SHEARDEFORMATION_SHELLSIMULATION_H

#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"
#include <Eigen/Dense>

class ShellSimulation
{
    
public:
    ShellSimulation() : _is_initialized(false), _itr_times(0) {};
    
    ~ShellSimulation() {};
    
public:
    void add_noise();
    
    bool set_up_simulation(const std::string &prefix, const std::string tar_shape);
    
    void energy_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df);
    /*
     This function calculate the overall engergy and its gradient at the given position x.
     @param[in] x:     The variable
     @param[in] ptr:   A void pointer, which is required by alglib
     
     @param[out] f:    The overall energy function
     @param[out] df:   The gradient of the overall energy function
     
     */
    
    void compute_deformed_surface(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
    /*
     This function is aimed to compute the deformed surface based on the interior energy and external energy
     defined in this class.
     f = E_S + E_external
     @param[in] V: The vertice of the surface
     @param[in] F: The triangle meshes
     
     @param[out] V: The vertice of the surface after deformation
     @param[out] F: The triangle meshes of the surface after deformation
     */
    
public:
    Eigen::MatrixXd VU;
    Eigen::MatrixXi F;
    double _PoissonsRatio;
    double _YoungsModulus;
    double _thickness;
    // Fixed points
    std::vector<Eigen::Vector3d> p_fixed;
    std::vector<int> p_fixed_index;
    std::vector<Eigen::Vector3d> external_force;
    double ratio;
    Eigen::VectorXd _omega_list;
    
private:
    bool _is_initialized;
    int _itr_times;
};


#endif //SHEARDEFORMATION_SHELLSIMULATION_H
