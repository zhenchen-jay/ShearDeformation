//
// Created by Zhen Chen on 12/27/18.
//
#pragma once

#include "ShellEnergy.h"

class ShellEnergyMatrixI : public ShellEnergy {

public:

    // Streching Energy

    virtual void streching_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F,
                          double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE) override;

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

    virtual void bending_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F,
                        double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE) override;

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

    virtual void test_bending_energy() override;


    

};


