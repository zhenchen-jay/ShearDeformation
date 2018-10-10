//
//  GravityEnergy.h
//  ShearDeformation
//
//  Created by Zhen Chen on 10/4/18.
//

#ifndef GravityEnergy_h
#define GravityEnergy_h
#include<Eigen/Dense>
#include<vector>

class GravityEnergy
{
public:
    GravityEnergy(){};
    ~GravityEnergy(){};
    
public:
    void gravity_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, std::vector<Eigen::Vector3d> external_force, double &E, Eigen::VectorXd &dE);
    
    
};


#endif /* GravityEnergy_h */
