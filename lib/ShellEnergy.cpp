//
// Created by Zhen Chen on 10/3/18.
//
#include <igl/vertex_triangle_adjacency.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <iomanip>
#include "ShellEnergy.h"
#include "GeometryFeature.h"

void ShellEnergy::streching_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F, double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE)
{
    double alpha = YoungsModulus*PossionRatio/((1+PossionRatio)*(1-2.0
            *PossionRatio));
    double beta = YoungsModulus/(2.0*(1+PossionRatio));
    E = 0;
    std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
    std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
    std::vector<Eigen::Matrix2d> IU_list; // The list of the first fundamental form of the original surface
    std::vector<Eigen::Matrix2d> ID_list;// The list of the first fundamental form of the current surface
    std::vector<double> dA_list;
    // get the information of vertex-triangle adjacency
    igl::vertex_triangle_adjacency(V0.rows(), F, VF, VFi);

    // Calculate the energy
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I_D, I_U;
        GeoFeature::calculate_first_fundamental_form(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),I_D);
        GeoFeature::calculate_first_fundamental_form(V0.row(F(i,0)),V0.row(F(i,1)),V0.row(F(i,2)),I_U);
        double dA = sqrt(I_U.determinant())/2.0; // 1/2 is the area of the parametric space
        dA_list.push_back(dA);
        IU_list.push_back(I_U);
        ID_list.push_back(I_D);

        Eigen::MatrixXd Q = I_U.inverse()*(I_D-I_U);
        E+=1.0/4.0 * thickness * dA * (alpha*0.5*(Q.trace()*Q.trace())+beta*(Q*Q).trace());
    }
    // Calculate the gradient
    dE.resize(3*V.rows());
    dE.setZero();
    for (int i=0;i<V.rows();i++)
    {

        for (int j = 0; j < VF[i].size(); j++) // Consider the faces which contains the vertex
        {
            Eigen::Matrix2d Q;
            std::vector<Eigen::Matrix2d> dI, dQ;
            int f = VF[i][j];
            int fi = VFi[i][j];
            GeoFeature::diff_first_fundamental_form(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), fi, dI);
            Q = IU_list[f].inverse() * (ID_list[f] - IU_list[f]);
            dQ.resize(3);
            for (int k = 0; k < 3; k++)
            {
                dQ[k] = IU_list[f].inverse() * dI[k];
                dE(3 * i + k) = dE(3 * i + k) + 1.0 / 2.0 * thickness * dA_list[f] *
                                                                (0.5 * alpha * Q.trace() * dQ[k].trace() +
                                                                 beta * (Q * dQ[k]).trace());
            }
        }
    }
    // std::cout<<dE.norm()<<std::endl;
}



void ShellEnergy::bending_energy(Eigen::MatrixXd V, Eigen::MatrixXd V0, Eigen::MatrixXi F, double YoungsModulus, double PossionRatio, double thickness, double &E, Eigen::VectorXd &dE)
{
    
}
