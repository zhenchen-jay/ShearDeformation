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
    if (start_ind != 0 && start_ind != 1 && start_ind != 2 && start_ind != -1)
    {
        std::cout << "The differential vertex should lie on the triangle or set start_ind as -1" << std::endl;
        return;
    }
    nf = ((V1-V0).cross(V2-V0)).normalized();
//    nf = ((V1-V0).cross(V2-V0));
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
    else if(start_ind==2)
    {
        e=V1-V0;
    }
    else
    {
        e << 0, 0 ,0;
    }
    for(int i=0;i<3;i++)
    {
        dn[i] = nf.dot(W.row(i)) * 1.0 / A * e.cross(nf);
//        Eigen::Vector3d normalized_nf = nf.normalized();
//        dn[i] = normalized_nf.dot(W.row(i)) * e.cross(normalized_nf);
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
    if (start_ind != 0 && start_ind != 1 && start_ind != 2 && start_ind != -1)
    {
        std::cout << "The differential vertex should lie on the triangle or set start_ind as -1" << std::endl;
        return;
    }
    Eigen::Vector3d A11, A12, A21, A22;
    if (start_ind == 0)
    {
        A11 = -2 * (V1 - V0);
        A12 = -(V1 + V2 - 2 * V0);
        A21 = -(V1 + V2 - 2 * V0);
        A22 = -2 * (V2 - V0);
    }
    else if (start_ind == 1)
    {
        A11 = 2 * (V1 - V0);
        A12 = V2 - V0;
        A21 = V2 - V0;
        A22 << 0, 0, 0;
    }
    else if (start_ind == 2)
    {
        A11 << 0, 0, 0;
        A12 = V1 - V0;
        A21 = V1 - V0;
        A22 = 2 * (V2 - V0);
    }
    else
    {
        A11 << 0, 0, 0;
        A12 << 0, 0, 0;
        A21 << 0, 0, 0;
        A22 << 0, 0, 0;
    }
    dI.resize(3);
    for(int k=0;k<3;k++)
    {
        dI[k] << A11(k), A12(k),
        A21(k), A22(k);
    }
    
}

void GeoFeature::calculate_second_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Vector3d V_0, Eigen::Vector3d V_1, Eigen::Vector3d V_2, std::vector<bool> real_pts, Eigen::Matrix2d &II)
{
    Eigen::Vector3d nf = ((V1-V0).cross(V2-V0));
    Eigen::Vector3d nf0, nf1, nf2;
    nf0.setZero();
    nf1.setZero();
    nf2.setZero();
    if(real_pts[0])
        nf0 = ((V_0-V1).cross(V2-V1)).normalized();
    if(real_pts[1])
        nf1 = ((V_1-V2).cross(V0-V2)).normalized();
    if(real_pts[2])
        nf2 = ((V_2-V0).cross(V1-V0)).normalized();
    
    Eigen::Vector3d n0 = (nf+nf0).normalized();
    Eigen::Vector3d n1 = (nf+nf1).normalized();
    Eigen::Vector3d n2 = (nf+nf2).normalized();
    
    II(0,0) = 2*(V1-V0).dot(n1-n0);
    II(0,1) = 2*(V1-V0).dot(n2-n0);
    II(1,0) = 2*(V2-V0).dot(n1-n0);
    II(1,1) = 2*(V2-V0).dot(n2-n0);
}


void GeoFeature::diff_second_fundamental_form(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Vector3d V_0, Eigen::Vector3d V_1, Eigen::Vector3d V_2, int start_ind, std::vector<bool> real_pts, std::vector<Eigen::Matrix2d> &dII)
{
    if(start_ind < 0 || start_ind > 6)
    {
        std::cout<<"The differential vertex should lie on the triangle or its adjacencies."<<std::endl;
        return;
    }
    if(start_ind<3) // Calculate the gradient w.r.t. the vertex lying on the triangle
    {
        Eigen::Matrix3d tri_vert; // start with V{start_ind}
        tri_vert.row((0-start_ind+3)%3) = V0;
        tri_vert.row((1-start_ind+3)%3) = V1;
        tri_vert.row((2-start_ind+3)%3) = V2;
        
        Eigen::Matrix3d adjacent_tri_vert;
        std::vector<bool> real_pts_after(3);
        adjacent_tri_vert.row((0-start_ind+3)%3) = V_0;
        adjacent_tri_vert.row((1-start_ind+3)%3) = V_1;
        adjacent_tri_vert.row((2-start_ind+3)%3) = V_2;
        
        real_pts_after[(0-start_ind+3)%3] = real_pts[0];
        real_pts_after[(1-start_ind+3)%3] = real_pts[1];
        real_pts_after[(2-start_ind+3)%3] = real_pts[2];
        
        Eigen::Vector3d nf;
        std::vector<Eigen::Vector3d> edge_n;
        std::vector<Eigen::Vector3d> dnf;
        std::vector<Eigen::Vector3d> adjacent_nf;
        std::vector<std::vector<Eigen::Vector3d> > diff_adjacent_nf(3, std::vector<Eigen::Vector3d>(3));
        std::vector<std::vector<Eigen::Vector3d> > diff_edge_n(3, std::vector<Eigen::Vector3d>(3));
        
        edge_n.resize(3);
        adjacent_nf.resize(3);
        
        face_normal(tri_vert.row(0), tri_vert.row(1), tri_vert.row(2), 0, nf, dnf);
        
        for(int i=0;i<3;i++)
        {
            if(real_pts_after[i])
                face_normal(tri_vert.row((i+2)%3), tri_vert.row((i+1)%3), adjacent_tri_vert.row(i), i-1 , adjacent_nf[i], diff_adjacent_nf[i]);
            else
            {
                adjacent_nf[i].setZero();
                diff_adjacent_nf[i].resize(3);
                for(int j=0;j<3;j++)
                    diff_adjacent_nf[i][j].setZero();
            }
        }
        
        for(int i=0;i<3;i++)
        {
            edge_n[i] = (nf+adjacent_nf[i]).normalized();
            for(int j=0;j<3;j++)
            {
                diff_edge_n[i][j] = 1.0 / (nf+adjacent_nf[i]).norm() * (dnf[j]+diff_adjacent_nf[i][j]-(dnf[j]+diff_adjacent_nf[i][j]).dot(edge_n[i])*edge_n[i]);
            }
        }
        
        std::vector<Eigen::Vector3d> e(3);
        e[0] << 1,0,0;
        e[1] << 0,1,0;
        e[2] << 0,0,1;
        dII.resize(3);
        
        for(int k=0;k<3;k++)
        {
            if(start_ind==0)
            {
                dII[k](0,0) = -2*e[k].dot(edge_n[1]-edge_n[0]) + 2*(tri_vert.row(1)-tri_vert.row(0)).dot(diff_edge_n[1][k]-diff_edge_n[0][k]);
                dII[k](0,1) = -2*e[k].dot(edge_n[2]-edge_n[0]) + 2*(tri_vert.row(1)-tri_vert.row(0)).dot(diff_edge_n[2][k]-diff_edge_n[0][k]);
                dII[k](1,0) = -2*e[k].dot(edge_n[1]-edge_n[0]) + 2*(tri_vert.row(2)-tri_vert.row(0)).dot(diff_edge_n[1][k]-diff_edge_n[0][k]);
                dII[k](1,1) = -2*e[k].dot(edge_n[2]-edge_n[0]) + 2*(tri_vert.row(2)-tri_vert.row(0)).dot(diff_edge_n[2][k]-diff_edge_n[0][k]);
            }
            else if(start_ind==1)
            {
                dII[k](0,0) = 2*e[k].dot(edge_n[0]-edge_n[2]) + 2*(tri_vert.row(0)-tri_vert.row(2)).dot(diff_edge_n[0][k]-diff_edge_n[2][k]);
                dII[k](0,1) = 2*e[k].dot(edge_n[1]-edge_n[2]) + 2*(tri_vert.row(0)-tri_vert.row(2)).dot(diff_edge_n[1][k]-diff_edge_n[2][k]);
                dII[k](1,0) = 2*(tri_vert.row(1)-tri_vert.row(2)).dot(diff_edge_n[0][k]-diff_edge_n[2][k]);
                dII[k](1,1) = 2*(tri_vert.row(1)-tri_vert.row(2)).dot(diff_edge_n[1][k]-diff_edge_n[2][k]);
            }
            else if(start_ind==2)
            {
                dII[k](0,0) = 2*(tri_vert.row(2)-tri_vert.row(1)).dot(diff_edge_n[2][k]-diff_edge_n[1][k]);
                dII[k](0,1) = 2*(tri_vert.row(2)-tri_vert.row(1)).dot(diff_edge_n[0][k]-diff_edge_n[1][k]);
                dII[k](1,0) = 2*e[k].dot(edge_n[2]-edge_n[1]) + 2*(tri_vert.row(0)-tri_vert.row(1)).dot(diff_edge_n[2][k]-diff_edge_n[1][k]);
                dII[k](1,1) = 2*e[k].dot(edge_n[0]-edge_n[1]) + 2*(tri_vert.row(0)-tri_vert.row(1)).dot(diff_edge_n[0][k]-diff_edge_n[1][k]);
            }
            
        }
    }
    
    else // Calculate the gradient w.r.t. the vertex lying on the adjacent triangle
    {
        start_ind = start_ind - 3;
        Eigen::Matrix3d tri_vert; // start with V_{start_ind}
        tri_vert.row((0-start_ind+3)%3) = V0;
        tri_vert.row((1-start_ind+3)%3) = V1;
        tri_vert.row((2-start_ind+3)%3) = V2;
        
        Eigen::Matrix3d adjacent_tri_vert;
        adjacent_tri_vert.row((0-start_ind+3)%3) = V_0;
        adjacent_tri_vert.row((1-start_ind+3)%3) = V_1;
        adjacent_tri_vert.row((2-start_ind+3)%3) = V_2;
        
        std::vector<bool> real_pts_after(3);
        real_pts_after[(0-start_ind+3)%3] = real_pts[0];
        real_pts_after[(1-start_ind+3)%3] = real_pts[1];
        real_pts_after[(2-start_ind+3)%3] = real_pts[2];
        
        Eigen::Vector3d nf;
        std::vector<Eigen::Vector3d> edge_n;
        std::vector<Eigen::Vector3d> dnf;
        std::vector<Eigen::Vector3d> adjacent_nf;
        std::vector<std::vector<Eigen::Vector3d> > diff_adjacent_nf(3, std::vector<Eigen::Vector3d>(3));
        std::vector<std::vector<Eigen::Vector3d> > diff_edge_n(3, std::vector<Eigen::Vector3d>(3));
        
        edge_n.resize(3);
        adjacent_nf.resize(3);
        
        face_normal(tri_vert.row(0), tri_vert.row(1), tri_vert.row(2), -1, nf, dnf);
        
        for(int i=0;i<3;i++)
        {
            if(real_pts_after[i])
                face_normal(tri_vert.row((i+2)%3), tri_vert.row((i+1)%3), adjacent_tri_vert.row(i), 2, adjacent_nf[i], diff_adjacent_nf[i]);
            else
            {
                adjacent_nf[i].setZero();
                diff_adjacent_nf[i].resize(3);
                for(int j=0;j<3;j++)
                    diff_adjacent_nf[i][j].setZero();
            }
        }
        
        for(int i=0;i<3;i++)
        {
            edge_n[i] = (nf+adjacent_nf[i]).normalized();
            for(int j=0;j<3;j++)
            {
                diff_edge_n[i][j] = 1.0 / (nf+adjacent_nf[i]).norm() * (dnf[j]+diff_adjacent_nf[i][j]-(dnf[j]+diff_adjacent_nf[i][j]).dot(edge_n[i])*edge_n[i]);
            }
        }
        
        dII.resize(3);
        
        for(int k=0;k<3;k++)
        {
            if(start_ind==0)
            {
                dII[k](0,0) = 2*(tri_vert.row(1)-tri_vert.row(0)).dot(-diff_edge_n[0][k]);
                dII[k](0,1) = 2*(tri_vert.row(1)-tri_vert.row(0)).dot(-diff_edge_n[0][k]);
                dII[k](1,0) = 2*(tri_vert.row(2)-tri_vert.row(0)).dot(-diff_edge_n[0][k]);
                dII[k](1,1) = 2*(tri_vert.row(2)-tri_vert.row(0)).dot(-diff_edge_n[0][k]);
            }
            else if(start_ind==1)
            {
                dII[k](0,0) = 2*(tri_vert.row(0)-tri_vert.row(2)).dot(diff_edge_n[0][k]);
                dII[k](0,1) = 0;
                dII[k](1,0) = 2*(tri_vert.row(1)-tri_vert.row(2)).dot(diff_edge_n[0][k]);
                dII[k](1,1) = 0;
            }
            else if(start_ind==2)
            {
                dII[k](0,0) = 0;
                dII[k](0,1) = 2*(tri_vert.row(2)-tri_vert.row(1)).dot(diff_edge_n[0][k]);
                dII[k](1,0) = 0;
                dII[k](1,1) = 2*(tri_vert.row(0)-tri_vert.row(1)).dot(diff_edge_n[0][k]);
            }
            
        }
    }
    
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


void GeoFeature::test_second_fundamental_form()
{
    Eigen::Matrix3d triangle, adjacent_tiangle,triangle1,adjacent_tiangle1;
    triangle<<1,1,0,
    -1,1,0,
    -1,-1,0;
    adjacent_tiangle<<0,0,0,
    1,-1,-0.5,
    0,0,-1.5;
    triangle1 = triangle;
    adjacent_tiangle1 = adjacent_tiangle;
    Eigen::Matrix2d II,II_1,diff_II;
    std::vector<Eigen::Matrix2d> dII;
    std::vector<bool> real_pts(3);
    for(int i=0;i<3;i++)
        real_pts[i] = false;
    real_pts[1] = true;
    
    calculate_second_fundamental_form(triangle.row(0), triangle.row(1), triangle.row(2), adjacent_tiangle.row(0), adjacent_tiangle.row(1), adjacent_tiangle.row(2), real_pts,II);
    
    std::cout<<"Test for the vertices lying on the current triabgle"<<std::endl;
    // Add some disturbance to the vertex of the triangle
    for(int i=0;i<3;i++)
    {
        diff_second_fundamental_form(triangle.row(0), triangle.row(1), triangle.row(2), adjacent_tiangle.row(0), adjacent_tiangle.row(1), adjacent_tiangle.row(2), i, real_pts,dII);
        Eigen::Vector3d w;
        w.setZero();
        for(int j=0;j<3;j++)
        {
            for(int k=3;k<9;k++)
            {
                w(j) = pow(10,-k);
                triangle1(i,j) = triangle(i,j) + w(j);
                calculate_second_fundamental_form(triangle1.row(0), triangle1.row(1), triangle1.row(2), adjacent_tiangle.row(0), adjacent_tiangle.row(1), adjacent_tiangle.row(2), real_pts,II_1);
                diff_II = (II_1-II)/w(j);
                std::cout<<"Eplison is: "<<w(j)<<std::endl;
                std::cout<<"The finite difference of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<diff_II<<std::endl;
                std::cout<<"The gradient of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<dII[j]<<std::endl;
                std::cout<<"The error between the gradient and the finite difference of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<(diff_II-dII[j]).cwiseAbs()<<std::endl;
            }
            w.setZero();
            triangle1 = triangle;
        }
    }
    
    std::cout<<"Test for the vertices lying on the adjacent triangles"<<std::endl;
    // Add some disturbance to the vertex of adjacent triangles
    for(int i=0;i<3;i++)
    {
        diff_second_fundamental_form(triangle.row(0), triangle.row(1), triangle.row(2), adjacent_tiangle.row(0), adjacent_tiangle.row(1), adjacent_tiangle.row(2), i+3, real_pts,dII);
        Eigen::Vector3d w;
        w.setZero();
        for(int j=0;j<3;j++)
        {
            for(int k=3;k<9;k++)
            {
                w(j) = pow(10,-k);
                adjacent_tiangle1(i,j) = adjacent_tiangle(i,j) + w(j);
                calculate_second_fundamental_form(triangle.row(0), triangle.row(1), triangle.row(2), adjacent_tiangle1.row(0), adjacent_tiangle1.row(1), adjacent_tiangle1.row(2), real_pts, II_1);
                diff_II = (II_1-II)/w(j);
                std::cout<<"Eplison is: "<<w(j)<<std::endl;
                std::cout<<"The finite difference of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<diff_II<<std::endl;
                std::cout<<"The gradient of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<dII[j]<<std::endl;
                std::cout<<"The error between the gradient and the finite difference of the II w.r.t. "<<j<<"th component of the vertex "<<i<<" is: "<<std::endl;
                std::cout<<(diff_II-dII[j]).cwiseAbs()<<std::endl;
            }
            w.setZero();
            adjacent_tiangle1 = adjacent_tiangle;
        }
    }
    
}
