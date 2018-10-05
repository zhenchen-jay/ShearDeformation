//
//  GeometryFeature.cpp
//  ShearDeformation_bin
//
//  Created by Zhen Chen on 10/5/18.
//

#include <iostream>
#include "GeometryFeature.h"


void GeoFeature::face_normal(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int start_ind, Eigen::Vector3d &nf, std::vector<Eigen::Vector3d> &dn)
{
    if (start_ind != 0 && start_ind != 1 && start_ind != 2)
    {
        std::cout << "The differential vertex should lie on the triangle" << std::endl;
        return;
    }
    nf = ((V1-V0).cross(V2-V0)).normalized();
    Eigen::Matrix3d W;
    W << 1,0,0,
    0,1,0,
    0,0,1;
    double A = (V1-V0).cross(V2-V0).norm();
    dn.resize(3);
    Eigen::Vector3d e;
    if(start_ind==0)
    {
        e=V2-V1;
    }
    else if(start_ind==1)
    {
        e=V0-V2;
    }
    else
    {
        e=V1-V0;
    }
    for(int i=0;i<3;i++)
    {
        dn[i] = nf.dot(W.row(i)) * 1.0 / A * e.cross(nf);
    }
    
}

void GeoFeature::calculate_first_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2,Eigen::Matrix2d &I)
{
    I(0,0) = (V1-V0).dot(V1-V0);
    I(1,0) = (V1-V0).dot(V2-V0);
    I(0,1) = (V1-V0).dot(V2-V0);
    I(1,1) = (V2-V0).dot(V2-V0);
}

void GeoFeature::diff_first_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int start_ind,std::vector<Eigen::Matrix2d> &dI)
{
    if (start_ind != 0 && start_ind != 1 && start_ind != 2)
    {
        std::cout << "The differential vertex should lie on the triangle" << std::endl;
        return;
    }
    Eigen::Vector3d A11, A12, A21, A22;
    if (start_ind == 0)
    {
        A11 = -2 * (V1 - V0);
        A12 = -(V1 + V2 - 2 * V0);
        A21 = -(V1 + V2 - 2 * V0);
        A22 = -2 * (V2 - V0);
    } else if (start_ind == 1)
    {
        A11 = 2 * (V1 - V0);
        A12 = V2 - V0;
        A21 = V2 - V0;
        A22 << 0, 0, 0;
    } else
    {
        A11 << 0, 0, 0;
        A12 = V1 - V0;
        A21 = V1 - V0;
        A22 = 2 * (V2 - V0);
    }
    dI.resize(3);
    for(int k=0;k<3;k++)
    {
        dI[k] << A11(k), A12(k),
        A21(k), A22(k);
    }
    
}

void GeoFeature::calculate_second_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Vector3d n0, Eigen::Vector3d n1, Eigen::Vector3d n2, Eigen::Matrix2d &II)
{
    II(0,0) = 2*(V1-V0).dot(n1-n0);
    II(0,1) = 2*(V1-V0).dot(n2-n0);
    II(1,0) = 2*(V2-V0).dot(n1-n0);
    II(1,1) = 2*(V2-V0).dot(n2-n0);
}








void GeoFeature::test_face_normal()
{
    Eigen::Matrix3d V;
    Eigen::Matrix3d U;
    Eigen::Vector3d n,n1;
    std::vector<Eigen::Vector3d> dn;
    std::vector<Eigen::Vector3d> dn1;
    
    V << 1,0,0,
    0,0.5,0,
    0,0,2;
    
    U=V;
    
    Eigen::Vector3d w;
    w.setZero();
    
    for(int j=0;j<3;j++)
    {
        face_normal(V.row(0), V.row(1), V.row(2), j, n, dn);
        for(int i=0;i<3;i++)
        {
            for(int k=3;k<7;k++)
            {
                w(i) = pow(10,-k);
                U(j,i) = V(j,i)+w(i);
                face_normal(U.row(0), U.row(1), U.row(2), j, n1, dn1);
                
                Eigen::Vector3d diff_n = (n1-n) / w(i);
                std::cout<<"The test vertex is V"<<j<<std::endl;
                std::cout<<"The variable is the "<<i<<"th component"<<std::endl;
                std::cout<<"Eplison is: "<<w(i)<<std::endl;
                Eigen::Vector3d grad = (dn[0]*w(0)+dn[1]*w(1)+dn[2]*w(2))/w(i);
                for(int p=0;p<3;p++)
                {
                    std::cout<<"The "<<p<<"th component of the finite differecce is: "<<diff_n(p)<<std::endl;
                    std::cout<<"The "<<p<<"th component of the gradient is: "<<grad(p)<<std::endl;
                    std::cout<<"The difference is: "<<abs(diff_n(p)-grad(p))<<std::endl;
                    
                }
            }
            w.setZero();
            U=V;
        }
    }
    
}
