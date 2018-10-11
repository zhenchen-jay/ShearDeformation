//
//  GravityEnergy.cpp
//  ShearDeformation_bin
//
//  Created by 陈臻 on 10/4/18.
//

#include "GravityEnergy.h"
#include <iostream>

void GravityEnergy::gravity_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, std::vector<Eigen::Vector3d> external_force, double &E, Eigen::VectorXd &dE)
{
    E = 0;
    dE.resize(3*V.rows());
    dE.setZero();
    if(external_force.size()!=V.rows())
    {
        std::cout<<"The number of external forces doesn't match the the mumber of the vertices"<<std::endl;
        return;
        
    }
    for(int i=0;i<V.rows();i++)
    {
        E -= external_force[i].dot(V.row(i)-V0.row(i));
        dE.segment(3*i, 3) = -external_force[i];
    }
}
