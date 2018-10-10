//
// Created by 陈臻 on 10/3/18.
//

#ifndef SHEARDEFORMATION_SHELLENERGY_H
#define SHEARDEFORMATION_SHELLENERGY_H

#endif //SHEARDEFORMATION_SHELLENERGY_H

#ifndef SHEARDEFORMATION_SIMULATION_H
#define SHEARDEFORMATION_SIMULATION_H
#include <Eigen/Dense>
#include <vector>
#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"

class ShellEnergy{
public:
    ShellEnergy() {};
    
    ~ShellEnergy() {};
    
public:
    
    // Streching Energy
    
    void streching_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F,
                          double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE);
    
    /*
     This function calculates the streching energy
     @param[in] V:  The vertices of the surface
     @param[in] V0: The vertices of the original surface
     @param[in] YoungsModulus: Young's Modulus
     @param[in] PossionRatio: Possion Ratio
     @param[in] Thickness: The thickness of the shell
     
     @parma[out] E: The streching Energy
     @parma[out] dE: The gradient of the streching Energy
     */
    
    // Bending Energy
    
    void bending_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F,
                        double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE);
    
    /*
     This function calculates the bending energy
     @param[in] V:  The vertices of the surface
     @param[in] V0: The vertices of the original surface
     @param[in] YoungsModulus: Young's Modulus
     @param[in] PossionRatio: Possion Ratio
     @param[in] Thickness: The thickness of the shell
     
     @parma[out] E: The bending Energy
     @parma[out] dE: The gradient of the bending Energy
     */
    
    void test_bending_energy();
};


#endif //SHEARDEFORMATION_SIMULATION_H
