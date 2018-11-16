//
// Created by Zhen Chen on 11/2/18.
//

#ifndef PROJECT_COMPUTECOEFFICIENT_H
#define PROJECT_COMPUTECOEFFICIENT_H
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"


class ComputeCoefficient
{
public:
    ComputeCoefficient(){}
    ~ComputeCoefficient(){}
    
public:
    Eigen::VectorXd find_optimal_sol(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, double YonungsModulus, double PoissonRatio, double thickness, double M1, double M2);
    
    Eigen::VectorXd get_coefficient(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, double YonungsModulus, double PoissonRatio, double thickness, double M1, double M2);
    /*
     This funciton is used to compute the coefficients of the first fundamental form.
     @param[in]: V0: The vertices of the undeformed surface
     @param[in]: V: The vertices of the target surface
     @param[in]: F: Geometric connection of the vertices, assuming that no flip happened
     @param[in]: {YoungsModulus, PoissonRatio, thickness}: The material parameters
     @param[in]: M1: The lower bund of the coefficients
     @param[in]: M2: The upper bounf of the coefficients
     
     @return:   The vector which contains the coefficient information of the first fundamental form.
     
     */
    
    void test();
    void test_libigl();
    
private:
    //void construct_LS(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, double YoungsModulus, double PoissonRatio, double thickness);
    void construct_LS(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &A, Eigen::VectorXd &b, double YoungsModulus, double PoissonRatio, double thickness);

};
#endif //PROJECT_COMPUTECOEFFICIENT_H
