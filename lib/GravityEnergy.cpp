//
//  GravityEnergy.cpp
//  ShearDeformation_bin
//
//  Created by 陈臻 on 10/4/18.
//
#ifndef ORIGIN
#define ORIGIN -0.5
#endif

#include "GravityEnergy.h"

void GravityEnergy::gravity_energy(Eigen::MatrixXd V, double mass, double &E, Eigen::VectorXd &dE)
{
    E = 0;
    dE.resize(V.rows()*3);
    dE.setZero();
    for(int i=0;i<V.rows();i++)
    {
        if(V(i,1)-ORIGIN>0)
        {
            E = E + 9.8*(V(i,1)-ORIGIN);
            dE(3*i+1) = 9.8;
        }
    }
}
