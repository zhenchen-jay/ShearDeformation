//
//  GravityEnergy.h
//  ShearDeformation
//
//  Created by Zhen Chen on 10/4/18.
//

#ifndef ExternalEnergy_h
#define ExternalEnergy_h
#include<Eigen/Dense>
#include<vector>

class ExternalEnergy
{
public:
    ExternalEnergy(){};
    ~ExternalEnergy(){};
    
public:
    void external_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, std::vector<Eigen::Vector3d> external_force, double &E, Eigen::VectorXd &dE);
    
    
};


#endif /* GravityEnergy_h */
