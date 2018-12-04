#include <unsupported/Eigen/MatrixFunctions>
#include <igl/readOBJ.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include "FindFirstFundamentalCoef.h"
#include "GeometryFeature.h"
#include "MeshConnection.h"
#include "ShellEnergyStandard.h"

void FindFirstFundamentalCoef::set_up(Eigen::MatrixXd VD, Eigen::MatrixXi F0, double YoungsModulus, double PoissonRatio, double thickness)
{
    V = VD;
    F = F0;
    _YoungsModulus = YoungsModulus;
    _PoissonsRatio = PoissonRatio;
    _thickness = thickness;
    _alpha = _YoungsModulus*_PoissonsRatio/((1+_PoissonsRatio)*(1-2.0
                                                                *_PoissonsRatio));
    _beta = _YoungsModulus/(2.0*(1+_PoissonsRatio));
    if(!_is_initialized)
    {
        std::cout<<"Start initialization"<<std::endl;
        MeshConnection::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
        MeshConnection::triangle_triangle_adjacency(F, TT, TTi);
        ID_list.assign(F.rows(), Eigen::Matrix2d());
        IID_list.assign(F.rows(), Eigen::Matrix2d());
        
        // Get the basic variables
        std::cout<<"Initialization: get the basic variables"<<std::endl;
        std::vector<int> index(3);
        std::vector<Eigen::Vector3d> adjacent_points(3);
        std::vector<Eigen::Vector3d> adjacent_points_U(3);
        std::vector<bool> real_pts(3);
        for(int i=0;i<F.rows();i++)
        {
            // Fisrt fundamental form
            Eigen::Matrix2d I_D, II_D;
            GeoFeature::calculate_first_fundamental_form(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),I_D);
            ID_list[i]=I_D;
            
            
            // Second fundamental form
            for(int k=0;k<3;k++)
            {
                if(TT(i,k)!=-1)
                {
                    index[(k+2)%3] = F( TT(i,k), (TTi(i,k)+2)%3 );
                    adjacent_points[(k+2)%3] = V.row(index[(k+2)%3]);
                    real_pts[(k+2)%3] = true;
                }
                else
                {
                    index[(k+2)%3] = -1; // Boundary edge
                    adjacent_points[(k+2)%3] << -1,-1,-1;
                    adjacent_points_U[(k+2)%3] << -1,-1,-1;
                    real_pts[(k+2)%3] = false;
                }
            }
            
            GeoFeature::calculate_second_fundamental_form(V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), adjacent_points[0], adjacent_points[1], adjacent_points[2], real_pts, II_D);
            IID_list[i]=II_D;
        }
        dID_list.assign(V.rows(), std::vector<Eigen::Matrix2d >{});
        // Get the gradient information of the first and second fundamental form
        std::cout<<"Initialization: get the gradient information of the first fundamental form"<<std::endl;
        for(int i=0;i<V.rows();i++)
        {
            for (int j = 0; j < VF[i].size(); j++) // Consider the faces which contains the vertex
            {
                std::vector<Eigen::Matrix2d> dI;
                int f = VF[i][j];
                int fi = VFi[i][j];
                GeoFeature::diff_first_fundamental_form(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), fi, dI);
                for(int r = 0;r<3;r++)
                {
                    dID_list[i].push_back(dI[r]);
                }
            }
        }
        std::cout<<"Initialization: getting the gradient information of the first fundamental form finished!"<<std::endl;
        dIID_list_neighbor.assign(V.rows(), std::vector<Eigen::Matrix2d >{});
        dIID_list_2_neighbor.assign(V.rows(), std::vector<Eigen::Matrix2d >{});
        std::cout<<"Initialization: get the gradient information of the second fundamental form"<<std::endl;
        for(int i=0;i<V.rows();i++)
        {
            for(int j=0;j<VF[i].size();j++) // Consider the faces which contains the vertex
            {
                std::vector<Eigen::Matrix2d> dII;
                int f = VF[i][j];
                int fi = VFi[i][j];
                for(int k=0;k<3;k++)
                {
                    if(TT(f,k)!=-1)
                    {
                        index[(k+2)%3] = F( TT(f,k), (TTi(f,k)+2)%3 );
                        adjacent_points[(k+2)%3] = V.row(index[(k+2)%3]);
                        real_pts[(k+2)%3] = true;
                    }
                    else
                    {
                        index[(k+2)%3] = -1; // Boundary edge
                        adjacent_points[(k+2)%3] << -1,-1,-1;
                        real_pts[(k+2)%3] = false;
                    }
                }
                GeoFeature::diff_second_fundamental_form(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), adjacent_points[0], adjacent_points[1], adjacent_points[2], fi, real_pts, dII);
                
                for(int r = 0;r<3;r++)
                {
                    dIID_list_neighbor[i].push_back(dII[r]);
                }
                
                if(TT(f,(fi+1)%3)!=-1) // Doesn't belong to the boundary
                {
                    int fa = TT(f,(fi+1)%3);
                    // We need to consider the gradient of II_{fa} w.r.t V
                    int start_ind = (TTi(f,(fi+1)%3) + 2) % 3 + 3;
                    for(int k=0;k<3;k++)
                    {
                        if(TT(fa,k)!=-1)
                        {
                            index[(k+2)%3] = F( TT(fa,k), (TTi(fa,k)+2)%3 );
                            adjacent_points[(k+2)%3] = V.row(index[(k+2)%3]);
                            real_pts[(k+2)%3] = true;
                        }
                        else
                        {
                            index[(k+2)%3] = -1; // Boundary edge
                            adjacent_points[(k+2)%3] << -1,-1,-1;
                            real_pts[(k+2)%3] = false;
                        }
                    }
                    GeoFeature::diff_second_fundamental_form(V.row(F(fa, 0)), V.row(F(fa, 1)), V.row(F(fa, 2)), adjacent_points[0], adjacent_points[1], adjacent_points[2], start_ind, real_pts, dII);
                    for(int r = 0;r<3;r++)
                    {
                        dIID_list_2_neighbor[i].push_back(dII[r]);
                    }
                }
            }
        }
        std::cout<<"Initialization: getting the gradient information of the second fundamental form finished!"<<std::endl;
    }
    _is_initialized = true;
    std::cout<<"Initialization finished!"<<std::endl;
}

void FindFirstFundamentalCoef::get_func_grad(Eigen::VectorXd &x, double &f, Eigen::VectorXd &df)
{
    // Valiadation check
    if(!_is_initialized)
    {
        std::cout<<"Please call set_up() function to set up first"<<std::endl;
        return;
    }
    if(x.size() != 3*F.rows())
    {
        std::cout<<"The size of the variables doesn't match F.rows()"<<std::endl;
        std::cout<<"x.size: "<<x.size()<<std::endl;
        std::cout<<"3*F.rows: "<<3*F.rows()<<std::endl;
        return;
    }
    
    f = 0;
    Eigen::VectorXd vec_f(3*V.rows()); // f = 1/2*vec_f.transpose()*vec_f
    vec_f.setZero();
    
    std::vector<Eigen::Matrix2d> IU_list; // The list of the first fundamental form of the original surface
    
    std::vector<double> dA_list;
    // get the information of vertex-triangle adjacency
    
    // Get the basic variables
    std::cout<<"Get the basic variables"<<std::endl;
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I_U;
        I_U << x(3*i),x(3*i+1),
        x(3*i+1),x(3*i+2);
        if(I_U.determinant()<=0)
        {
            std::cout<<"The first fundamental form isn't positive"<<std::endl;
            return;
        }
        double dA = sqrt(I_U.determinant())/2.0; // 1/2 is the area of the parametric space
        dA_list.push_back(dA);
        IU_list.push_back(I_U);
    }
    // Get the funciton
    std::cout<<"Get the funciton"<<std::endl;
    
    for (int i=0;i<V.rows();i++)
    {
        for (int j = 0; j < VF[i].size(); j++) // Consider the faces which contains the vertex
        {
            Eigen::Matrix2d Q;
            std::vector<Eigen::Matrix2d> dQ, dI;
            int f = VF[i][j];
            int fi = VFi[i][j];
            Q = IU_list[f].inverse() * (ID_list[f] - IU_list[f]);
            dQ.resize(3);
            GeoFeature::diff_first_fundamental_form(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), fi, dI);
            for (int k = 0; k < 3; k++)
            {
                dQ[k] = IU_list[f].inverse() * dID_list[i][3*j+k];
                vec_f(3 * i + k) = vec_f(3 * i + k) + 1.0 / 2.0 * _thickness * dA_list[f] *
                (0.5 * _alpha * Q.trace() * dQ[k].trace() +
                 _beta * (Q * dQ[k]).trace());
            }
        }
    }
    
    std::cout<<"First Fundamental form finished!"<<std::endl;
    std::vector<int> index(3);
    std::vector<Eigen::Vector3d> adjacent_points(3);
    std::vector<Eigen::Vector3d> adjacent_points_U(3);
    std::vector<bool> real_pts(3);
    //
    std::cout<<"Second fundamental form started"<<std::endl;
    for(int i=0;i<V.rows();i++)
    {
        int adj_cout = 0;
        for(int j=0;j<VF[i].size();j++) // Consider the faces which contains the vertex
        {
            Eigen::Matrix2d Q;
            std::vector<Eigen::Matrix2d> dQ;
            int f = VF[i][j];
            int fi = VFi[i][j];
            Q = IU_list[f].inverse()*IID_list[f];
            dQ.resize(3);
            for(int k=0;k<3;k++)
            {
                dQ[k] = IU_list[f].inverse()*dIID_list_neighbor[i][3*j+k];
                vec_f(3 * i + k) = vec_f(3 * i + k) + 1.0 / 6.0 * pow(_thickness,3.0) * dA_list[f] *
                (0.5 * _alpha * Q.trace() * dQ[k].trace() +
                 _beta * (Q * dQ[k]).trace());
                
            }
            
            if(TT(f,(fi+1)%3)!=-1) // Doesn't belong to the boundary
            {
                int fa = TT(f,(fi+1)%3);
                Q = IU_list[fa].inverse()*IID_list[fa];
                dQ.resize(3);
                for(int k=0;k<3;k++)
                {
                    dQ[k] = IU_list[fa].inverse()*dIID_list_2_neighbor[i][3*adj_cout+k];
                    vec_f(3 * i + k) = vec_f(3 * i + k) + 1.0 / 6.0 * pow(_thickness,3.0) * dA_list[fa] *
                    (0.5 * _alpha * Q.trace() * dQ[k].trace() +
                     _beta * (Q * dQ[k]).trace());
                }
                adj_cout ++;
            }
        }
    }
    std::cout<<"Second fundamental form finished"<<std::endl;
    f = 1.0/2*vec_f.transpose()*vec_f;
    
    // Get the gradient
    std::cout<<"Get the gradient"<<std::endl;
    Eigen::SparseMatrix<double> grad_vec_f(vec_f.size(),3*F.rows());
    grad_vec_f.setZero();
    Eigen::Matrix2d Id;
    Id.setIdentity();
    for(int i=0;i<V.rows();i++)
    {
        // The terms coming from the first fundamental form
        for(int j=0; j<VF[i].size();j++)
        {
            int f = VF[i][j];
            std::vector<Eigen::Matrix2d> dIU_inv(3),A(3),dI(3);
            std::vector<double> diff_sqrt_det(3);
            
            compute_derivative_sqrt_det(IU_list[f], diff_sqrt_det);
            compute_derivative_inv_mat(IU_list[f], dIU_inv);
            
            for(int r=0;r<3;r++)
            {
                for(int c=0;c<3;c++)
                {
                    double result = 0;
                    dI[r] << dID_list[i][3*j+r];
                    result += 1.0*_thickness/4.0*_alpha*( (dIU_inv[c]*ID_list[f]).trace() * (IU_list[f].inverse()*dI[r]).trace() + (IU_list[f].inverse()*ID_list[f]-Id).trace() * (dIU_inv[c]*dI[r]).trace() ) * sqrt(IU_list[f].determinant())*0.5;
                    result += 1.0*_thickness/4.0*_alpha*(IU_list[f].inverse()*ID_list[f]-Id).trace() * (IU_list[f].inverse()*dI[r]).trace()*diff_sqrt_det[c]*0.5;
                    result += 0.5*_thickness*_beta*(dIU_inv[c]*ID_list[f]*IU_list[f].inverse()*dI[r] + (IU_list[f].inverse()*ID_list[f]-Id)*dIU_inv[c]*dI[r]).trace()*sqrt(IU_list[f].determinant())*0.5;
                    result += 0.5*_thickness*_beta*( (IU_list[f].inverse()*ID_list[f]-Id)*IU_list[f].inverse()*dI[r] ).trace()*diff_sqrt_det[c]*0.5;
                    grad_vec_f.coeffRef(3*i+r,3*f+c) += result;
                }
            }
        }
        //        Eigen::MatrixXd grad_f_dMat;
        //        grad_f_dMat = Eigen::MatrixXd(grad_vec_f);
        //
        //        std::cout<<grad_f_dMat<<std::endl;
        int adj_count = 0;
        //        // The terms coming from the second fundamental form
        for(int j=0;j<VF[i].size();j++)
        {
            int f = VF[i][j];
            int fi = VFi[i][j];
            std::vector<Eigen::Matrix2d> dIU_inv(3),A(3),dII(3);
            std::vector<double> diff_sqrt_det(3);
            
            compute_derivative_sqrt_det(IU_list[f], diff_sqrt_det);
            compute_derivative_inv_mat(IU_list[f], dIU_inv);
            
            for(int r=0;r<3;r++)
            {
                for(int c=0;c<3;c++)
                {
                    double result = 0;
                    dII[r] << dIID_list_neighbor[i][3*j+r];
                    result += 1.0*pow(_thickness,3.0)/12.0*_alpha*( (dIU_inv[c]*IID_list[f]).trace() * (IU_list[f].inverse()*dII[r]).trace() + (IU_list[f].inverse()*IID_list[f]).trace() * (dIU_inv[c]*dII[r]).trace() ) * sqrt(IU_list[f].determinant())*0.5;
                    result += 1.0*pow(_thickness,3.0)/12.0*_alpha*(IU_list[f].inverse()*IID_list[f]).trace() * (IU_list[f].inverse()*dII[r]).trace()*diff_sqrt_det[c]*0.5;
                    result += 1.0*pow(_thickness,3.0)/6.0*_beta*(dIU_inv[c]*IID_list[f]*IU_list[f].inverse()*dII[r] + (IU_list[f].inverse()*IID_list[f])*dIU_inv[c]*dII[r]).trace()*sqrt(IU_list[f].determinant())*0.5;
                    result += 1.0*pow(_thickness,3.0)/6.0*_beta*( (IU_list[f].inverse()*IID_list[f])*IU_list[f].inverse()*dII[r] ).trace()*diff_sqrt_det[c]*0.5;
                    grad_vec_f.coeffRef(3*i+r,3*f+c) += result;
                }
            }
            
            if(TT(f,(fi+1)%3)!=-1)
            {
                int fa = TT(f,(fi+1)%3);
                compute_derivative_sqrt_det(IU_list[fa], diff_sqrt_det);
                compute_derivative_inv_mat(IU_list[fa], dIU_inv);
                for(int r=0;r<3;r++)
                {
                    for(int c=0;c<3;c++)
                    {
                        double result = 0;
                        dII[r] << dIID_list_2_neighbor[i][3*adj_count+r];
                        result += 1.0*pow(_thickness,3.0)/12.0*_alpha*( (dIU_inv[c]*IID_list[fa]).trace() * (IU_list[fa].inverse()*dII[r]).trace() + (IU_list[fa].inverse()*IID_list[fa]).trace() * (dIU_inv[c]*dII[r]).trace() ) * sqrt(IU_list[fa].determinant())*0.5;
                        result += 1.0*pow(_thickness,3.0)/12.0*_alpha*(IU_list[fa].inverse()*IID_list[fa]).trace() * (IU_list[fa].inverse()*dII[r]).trace()*diff_sqrt_det[c]*0.5;
                        result += 1.0*pow(_thickness,3.0)/6.0*_beta*(dIU_inv[c]*IID_list[fa]*IU_list[fa].inverse()*dII[r] + (IU_list[fa].inverse()*IID_list[fa])*dIU_inv[c]*dII[r]).trace()*sqrt(IU_list[fa].determinant())*0.5;
                        result += 1.0*pow(_thickness,3.0)/6.0*_beta*( (IU_list[fa].inverse()*IID_list[fa])*IU_list[fa].inverse()*dII[r] ).trace()*diff_sqrt_det[c]*0.5;
                        grad_vec_f.coeffRef(3*i+r,3*fa+c) += result;
                    }
                }
                adj_count++;
            }
        }
        
    }
    df = grad_vec_f.transpose()*vec_f;
    std::cout<<"Finished!"<<std::endl;
}

void FindFirstFundamentalCoef::compute_derivative_inv_mat(Eigen::Matrix2d A, std::vector<Eigen::Matrix2d> &dA)
{
    dA.resize(3);
    std::vector<Eigen::Matrix2d> C(3);
    C[0]<<0,0,
    0,1;
    C[1]<<0,-1,
    -1,0;
    C[2]<<1,0,
    0,0;
    
    dA[0]=-1.0*A(1,1)/A.determinant()*A.inverse() + 1.0/A.determinant()*C[0];
    dA[1]=2.0*A(0,1)/A.determinant()*A.inverse() + 1.0/A.determinant()*C[1];
    dA[2]=-1.0*A(0,0)/A.determinant()*A.inverse() + 1.0/A.determinant()*C[2];
    
    
}

void FindFirstFundamentalCoef::compute_derivative_sqrt_det(Eigen::Matrix2d A, std::vector<double> &diff_sqrt_det)
{
    diff_sqrt_det.resize(3);
    diff_sqrt_det[0] = 0.5*A(1,1)/sqrt(A.determinant());
    diff_sqrt_det[1] = -1.0*A(0,1)/sqrt(A.determinant());
    diff_sqrt_det[2] = 0.5*A(0,0)/sqrt(A.determinant());
}

void FindFirstFundamentalCoef::test_func_grad()
{
    double YoungsModulus = 1e5;
    double PossionRatio = 0.3;
    double thickness = 1e-3;
    
    Eigen::MatrixXd V0, V;
    Eigen::MatrixXi F0, F;
    
    igl::readOBJ("../../benchmarks/TestModels/rect.obj", V0, F0);
    igl::readOBJ("../../benchmarks/TestModels/saddle.obj", V, F);
    
//    igl::readOBJ("../../benchmarks/TestModels/test_squre.obj", V0,F0);
//    igl::readOBJ("../../benchmarks/TestModels/test_squre_bended.obj", V,F);
    double E, E1,E2,E3;
    Eigen::VectorXd dE, dE1,dE2,dE3;
    Eigen::VectorXd x(F.rows()*3);
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I;
        GeoFeature::calculate_first_fundamental_form(V0.row(F0(i,0)),V0.row(F0(i,1)),V0.row(F0(i,2)),I);
        x(3*i) = I(0,0);
        x(3*i+1) = I(0,1);
        x(3*i+2) = I(1,1);
    }
    set_up(V, F, YoungsModulus, PossionRatio, thickness);
    Eigen::VectorXd x1= x;
    get_func_grad(x, E, dE);
    auto op = std::make_unique<ShellEnergyStandard>();
    op->streching_energy(V, V0, F, YoungsModulus, PossionRatio, thickness, E2, dE2);
    op->bending_energy(V, V0, F, YoungsModulus, PossionRatio, thickness, E3, dE3);
    dE2=dE2 + dE3;
    E2 = E2 + E3;
    std::cout<<0.5*dE2.transpose()*dE2<<" "<<E<<std::endl;
    srand((unsigned)time(NULL));
    int selected_i = rand()%(F.rows()*3);
    Eigen::VectorXd eps(F.rows()*3);
    eps.setZero();
    for(int k=6;k<16;k++)
    {
        selected_i = 0;
        eps(selected_i) = pow(10,-k);
        x1 += eps;
        get_func_grad(x1, E1, dE1);
        std::cout<<"Selected index is: "<<selected_i<<" Eplison is: "<<eps(selected_i)<<std::endl;
        std::cout<<"finite difference is: "<<(E1-E)/pow(10,-k)<<std::endl;
        std::cout<<"gradient projection: "<<dE.dot(eps)/eps(selected_i)<<std::endl;
        std::cout<<"The different between the finite difference and gradient is: "<<abs((E1-E)/eps(selected_i) - dE.dot(eps)/eps(selected_i))<<std::endl;
        x1=x;
    }
    
    //    Eigen::Matrix2d test_A;
    //    std::vector<Eigen::Matrix2d> dA_list;
    //    test_A << 1.0, 0.5,
    //    0.5,2.0;
    //    compute_derative_inv_mat(test_A, dA_list);
    
}
