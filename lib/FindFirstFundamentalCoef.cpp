#include <unsupported/Eigen/MatrixFunctions>
#include <igl/readOBJ.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <memory>
#include <vector>
#include "FindFirstFundamentalCoef.h"
#include "GeometryFeature.h"
#include "MeshConnection.h"
#include "ShellEnergyStandard.h"

// TO DO: Calculate the second derivative and use newton method to solve this problem

void callback_func(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df, void* ptr)
{
    FindFirstFundamentalCoef* callback = (FindFirstFundamentalCoef*) ptr;
    callback->get_func_grad(x,f,df);
}

void FindFirstFundamentalCoef::compute_first_fundamental_form(Eigen::MatrixXd VD, Eigen::MatrixXi F0, std::vector<Eigen::Matrix2d> &IU_array, double YoungsModulus, double PoissonRatio, double thickness)
{
    _itr_times = 0;
    if(!_is_initialized)
    {
        set_up(VD, F0, YoungsModulus, PoissonRatio, thickness);
    }
    alglib::real_1d_array x;
    double epsg = 0;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0;
    alglib::ae_int_t maxits = 10000;
    //    alglib::mincgstate state;
    //    alglib::mincgreport rep;
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    x.setlength(3*F.rows());
    
    Eigen::MatrixXd V_int;
    Eigen::MatrixXi F_int;
    
    igl::readOBJ("../benchmarks/TestModels/rect.obj", V_int, F_int);
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I0;
        V_int = 0.5*V_int + 0.5*V;
        GeoFeature::calculate_first_fundamental_form(V_int.row(F_int(i,0)),V_int.row(F_int(i,1)),V_int.row(F_int(i,2)),I0);
        x[3*i] = sqrt(I0(0,0));
        x[3*i+1] = I0(0,1)/x[3*i];
        x[3*i+2] = sqrt(I0.determinant())/x[3*i];
    }
    
    // Solve the optimization problem using alglib
    //    alglib::mincgcreate(x, state);
    //    alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
    //    alglib::mincgsetstpmax(state, stpmax);
    //    alglib::mincgoptimize(state, callback_func, NULL, this);
    //    alglib::mincgresults(state, x, rep);
    alglib::minlbfgscreate(x.length(),x, state);
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    alglib::minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, callback_func, NULL, this);
    alglib::minlbfgsresults(state, x, rep);

    printf("%d\n", int(rep.terminationtype));
    std::ofstream outfile("L_list.dat",std::ios::trunc);
    for(int i=0;i<x.length()-1;i++)
    {
        outfile<<std::setprecision(16)<<x[i]<<"\n";
    }
    outfile<<std::setprecision(16)<<x[x.length()-1];
    outfile.close();
}

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
        //std::cout<<"Start initialization"<<std::endl;
        MeshConnection::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
        MeshConnection::triangle_triangle_adjacency(F, TT, TTi);
        ID_list.assign(F.rows(), Eigen::Matrix2d());
        IID_list.assign(F.rows(), Eigen::Matrix2d());
        
        // Get the basic variables
        //std::cout<<"Initialization: get the basic variables"<<std::endl;
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
        //std::cout<<"Initialization: get the gradient information of the first fundamental form"<<std::endl;
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
        //std::cout<<"Initialization: getting the gradient information of the first fundamental form finished!"<<std::endl;
        dIID_list_neighbor.assign(V.rows(), std::vector<Eigen::Matrix2d >{});
        dIID_list_2_neighbor.assign(V.rows(), std::vector<Eigen::Matrix2d >{});
        //std::cout<<"Initialization: get the gradient information of the second fundamental form"<<std::endl;
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
        //std::cout<<"Initialization: getting the gradient information of the second fundamental form finished!"<<std::endl;
    }
    _is_initialized = true;
    //std::cout<<"Initialization finished!"<<std::endl;
}

void FindFirstFundamentalCoef::get_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df)
{
    _itr_times++;
    if(_itr_times%100 == 0)
    {
        std::ofstream outfile("L_list.dat",std::ios::trunc);
        for(int i=0;i<x.length()-1;i++)
        {
            outfile<<std::setprecision(16)<<x[i]<<"\n";
        }
        outfile<<std::setprecision(16)<<x[x.length()-1];
        outfile.close();
    }
    // Valiadation check
    if(!_is_initialized)
    {
        std::cout<<"Please call set_up() function to set up first"<<std::endl;
        return;
    }
    if(x.length() != 3*F.rows())
    {
        std::cout<<"The size of the variables doesn't match F.rows()"<<std::endl;
        std::cout<<"x.size: "<<x.length()<<std::endl;
        std::cout<<"3*F.rows: "<<3*F.rows()<<std::endl;
        return;
    }
    
    f = 0;
    Eigen::VectorXd vec_f(3*V.rows()); // f = 1/2*vec_f.transpose()*vec_f
    vec_f.setZero();
    
    std::vector<Eigen::Matrix2d> IU_list; // The list of the first fundamental form of the original surface
    std::vector<Eigen::Matrix2d> L_list; // IU = L*L^T
    
    std::vector<double> dA_list;
    // get the information of vertex-triangle adjacency
    
    // Get the basic variables
    //std::cout<<"Get the basic variables"<<std::endl;
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I_U;
        Eigen::Matrix2d L;
        L << x[3*i],0,
        x[3*i+1],x[3*i+2];
        I_U = L*L.transpose();
        double dA = sqrt(I_U.determinant())/2.0; // 1/2 is the area of the parametric space
        dA_list.push_back(dA);
        IU_list.push_back(I_U);
        L_list.push_back(L);
    }
    // Get the funciton
    //std::cout<<"Get the funciton"<<std::endl;
    
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
    
    //std::cout<<"First Fundamental form finished!"<<std::endl;
    std::vector<int> index(3);
    std::vector<Eigen::Vector3d> adjacent_points(3);
    std::vector<Eigen::Vector3d> adjacent_points_U(3);
    std::vector<bool> real_pts(3);
    //
    //std::cout<<"Second fundamental form started"<<std::endl;
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
    //std::cout<<"Second fundamental form finished"<<std::endl;
    f = 1.0/2*vec_f.transpose()*vec_f;
    
    // Get the gradient
    //std::cout<<"Get the gradient"<<std::endl;
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

            compute_derivative_sqrt_det(L_list[f], diff_sqrt_det);
            compute_derivative_inv_mat(L_list[f], dIU_inv);

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
        // The terms coming from the second fundamental form
        for(int j=0;j<VF[i].size();j++)
        {
            int f = VF[i][j];
            int fi = VFi[i][j];
            std::vector<Eigen::Matrix2d> dIU_inv(3),A(3),dII(3);
            std::vector<double> diff_sqrt_det(3);

            compute_derivative_sqrt_det(L_list[f], diff_sqrt_det);
            compute_derivative_inv_mat(L_list[f], dIU_inv);

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
                compute_derivative_sqrt_det(L_list[fa], diff_sqrt_det);
                compute_derivative_inv_mat(L_list[fa], dIU_inv);
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
    Eigen::VectorXd gradf = grad_vec_f.transpose()*vec_f;
    for(int i=0;i<df.length();i++)
    {
        df[i] = gradf(i);
    }
    std::cout<<"Iteration Times: "<<_itr_times<<"\t"<<"Objective function: "<<f<<"\t"<<"||Gradient||: "<<gradf.norm()<<std::endl;
}

void FindFirstFundamentalCoef::compute_derivative_inv_mat(Eigen::Matrix2d A, std::vector<Eigen::Matrix2d> &dA)
{
    dA.resize(3);
    double x,y,z;
    x = A(0,0);
    y = A(1,0);
    z = A(1,1);
    Eigen::Matrix2d M;
    M << y*y+z*z,-x*y,
    -x*y,x*x;
    std::vector<Eigen::Matrix2d> C(3);
    C[0]<<0,-y,
    -y,2*x;
    C[1]<<2*y,-x,
    -x,0;
    C[2]<<2*z,0,
    0,0;
    
    dA[0]=1/(x*x*z*z)*C[0] - 2/(x*x*x*z*z)*M;
    dA[1]=1/(x*x*z*z)*C[1];
    dA[2]=1/(x*x*z*z)*C[2] - 2/(x*x*z*z*z)*M;
    
    
}

void FindFirstFundamentalCoef::compute_derivative_sqrt_det(Eigen::Matrix2d A, std::vector<double> &diff_sqrt_det)
{
    diff_sqrt_det.resize(3);
    int sign = 1;
    if(A.determinant()<0)
        sign = -1;
    diff_sqrt_det[0] = sign*A(1,1);
    diff_sqrt_det[1] = 0;
    diff_sqrt_det[2] = sign*A(0,0);
}

void FindFirstFundamentalCoef::test_func_grad()
{
    double YoungsModulus = 1e5;
    double PossionRatio = 0.3;
    double thickness = 1e-4;
    
    Eigen::MatrixXd V0, V;
    Eigen::MatrixXi F0, F;
    
    igl::readOBJ("../../benchmarks/TestModels/rect.obj", V0, F0);
    igl::readOBJ("../../benchmarks/TestModels/sphere.obj", V, F);
    
//    igl::readOBJ("../../benchmarks/TestModels/test_square.obj", V0,F0);
//    igl::readOBJ("../../benchmarks/TestModels/test_square_bended.obj", V,F);
    double E, E1,E2,E3;
    Eigen::VectorXd dE(F.rows()*3),dE2,dE3;
    alglib::real_1d_array alg_x,alg_x1,alg_df, alg_df1;
    Eigen::VectorXd x(F.rows()*3);
    alg_x.setlength(F.rows()*3);
    alg_df.setlength(F.rows()*3);
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I0;
        GeoFeature::calculate_first_fundamental_form(V0.row(F0(i,0)),V0.row(F0(i,1)),V0.row(F0(i,2)),I0);

        alg_x[3*i] = sqrt(I0(0,0));
        alg_x[3*i+1] = I0(0,1)/alg_x[3*i];
        alg_x[3*i+2] = sqrt(I0.determinant())/alg_x[3*i];
//        alg_x[3*i] = 1;
//        alg_x[3*i+1] = 0.03;
//        alg_x[3*i+2] = 0.05;
//        Eigen::Matrix2d L;
//        L << alg_x[3*i],0,
//        alg_x[3*i+1],alg_x[3*i+2];
//        std::cout<<L<<std::endl;
    }
    set_up(V, F, YoungsModulus, PossionRatio, thickness);
    Eigen::VectorXd x1= x;
    get_func_grad(alg_x, E, alg_df);
    for(int i=0;i<alg_df.length();i++)
    {
        dE(i) = alg_df[i];
    }
    auto op = std::make_unique<ShellEnergyStandard>();
    op->streching_energy(V, V0, F, YoungsModulus, PossionRatio, thickness, E2, dE2);
    op->bending_energy(V, V0, F, YoungsModulus, PossionRatio, thickness, E3, dE3);
    dE2=dE2 + dE3;
    E2 = E2 + E3;
    //std::cout<<0.5*dE2.transpose()*dE2<<" "<<E<<std::endl;
    srand((unsigned)time(NULL));
    int selected_i = rand()%(F.rows()*3);
    Eigen::VectorXd eps(F.rows()*3);
    eps.setZero();
    alg_x1 = alg_x;
    selected_i = 0;
    for(int k=3;k<11;k++)
    {
        eps(selected_i) = pow(10,-k);
        alg_x1[selected_i] += pow(10,-k);
        get_func_grad(alg_x1, E1, alg_df1);
        std::cout<<"Selected index is: "<<selected_i<<" Eplison is: "<<eps(selected_i)<<std::endl;
        std::cout<<"finite difference is: "<<(E1-E)/(pow(10,-k))<<std::endl;
        std::cout<<"gradient projection: "<<dE.dot(eps)/eps(selected_i)<<std::endl;
        std::cout<<"The different between the finite difference and gradient is: "<<std::abs((E1-E)/pow(10,-k) - dE.dot(eps)/eps(selected_i))<<std::endl;
        alg_x1 = alg_x;

    }
    
    //    Eigen::Matrix2d test_A;
    //    std::vector<Eigen::Matrix2d> dA_list;
    //    test_A << 1.0, 0.5,
    //    0.5,2.0;
    //    compute_derative_inv_mat(test_A, dA_list);
    
}
