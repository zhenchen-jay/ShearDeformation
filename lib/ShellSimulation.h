//
// Created by 陈臻 on 10/3/18.
//

#ifndef SHEARDEFORMATION_SHELLSIMULATION_H
#define SHEARDEFORMATION_SHELLSIMULATION_H

#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"
#include <Eigen/Dense>

class ShellSimulation
{

public:
    ShellSimulation() : _is_initialized(false) {};

    ~ShellSimulation() {};

public:

    void energy_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df);
    /*
    This function calculate the overall engergy and its gradient at the given position x.
    @param[in] x:     The variable
    @param[in] ptr:   A void pointer, which is required by alglib

    @param[out] f:    The overall energy function
    @param[out] df:   The gradient of the overall energy function

    */

    void compute_deformed_surface(Eigen::MatrixXd &V0, Eigen::MatrixXi &F0, double PossionRatio, double YoungsModulus, double thickness, double ratio);
    /*
    This function is aimed to compute the deformed surface based on the interior energy and external energy
    defined in this class.
    f = E_S + E_external
    @param[in] V0: The vertice of the surface
    @param[in] F0: The triangle meshes
    @oaram[in] PossionRatio: The Passion ratio
    @param[in] YoungsModulus: Young's Modulus
    @param[in] thickness:  Thickness
    @param[in] ratio: Streching ratio

    @param[out] V0: The vertice of the surface after deformation
    @param[out] F0: The triangle meshes of the surface after deformation
    */

private:
    Eigen::MatrixXd VU;
    Eigen::MatrixXi F;
    double _PossionRation;
    double _YoungsModulus;
    double _thickness;
    bool _is_initialized;
    // Fixed points
    std::vector<Eigen::Vector3d> p_fixed;
    std::vector<int> p_fixed_index;
};


#endif //SHEARDEFORMATION_SHELLSIMULATION_H
