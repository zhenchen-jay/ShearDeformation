//
//  GravityEnergy.h
//  ShearDeformation
//
//  Created by Zhen Chen on 10/4/18.
//

#ifndef GravityEnergy_h
#define GravityEnergy_h
#include<Eigen/Dense>

class GravityEnergy
{
public:
    GravityEnergy(){};
    ~GravityEnergy(){};
    
public:
    void gravity_energy(Eigen::MatrixXd V, double mass, double &E, Eigen::VectorXd &dE);
    
    
};


#endif /* GravityEnergy_h */
