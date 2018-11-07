//
// Created by Zhen Chen on 11/2/18.
//
#include <igl/readOBJ.h>
#include <iostream>
#include <cmath>
#include "ComputeCoefficient.h"
#include "GeometryFeature.h"
#include "MeshConnection.h"
#include "eiquadprog.h"

void ComputeCoefficient::construct_LS(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &A, Eigen::VectorXd &b, double YoungsModulus, double PoissonRatio, double thickness)
{
    double alpha = YoungsModulus*PoissonRatio/((1+PoissonRatio)*(1-2.0
                                                                 *PoissonRatio));
    double beta = YoungsModulus/(2.0*(1+PoissonRatio));
    std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
    std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
    Eigen::MatrixXi TT;             // The matrix which stores the information of adjacent faces
    Eigen::MatrixXi TTi;
    std::vector<Eigen::Matrix2d> IU_list; // The list of the first fundamental form of the original surface
    std::vector<Eigen::Matrix2d> ID_list;// The list of the first fundamental form of the current surface
    // std::vector<Eigen::Matrix2d> IIU_list; // The list of the first fundamental form of the current surface
    std::vector<Eigen::Matrix2d> IID_list;// The list of the second fundamental form of the current surface
    
    std::vector<int> index(3);
    std::vector<Eigen::Vector3d> adjacent_points(3);
    std::vector<Eigen::Vector3d> adjacent_points_U(3);
    std::vector<bool> real_pts(3);
    
    std::vector<double> dA_list;
    // get the information of vertex-triangle adjacency
    MeshConnection::vertex_triangle_adjacency(V0.rows(), F, VF, VFi);
    MeshConnection::triangle_triangle_adjacency(F, TT, TTi);
    
    A.resize(3*V.rows(), F.rows());
    b.resize(3*V.rows());
    A.setZero();
    b.setZero();
    
    // Calculate the energy
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I_D, I_U, II_D, II_U;
        GeoFeature::calculate_first_fundamental_form(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),I_D);
        GeoFeature::calculate_first_fundamental_form(V0.row(F(i,0)),V0.row(F(i,1)),V0.row(F(i,2)),I_U);
        double dA = sqrt(I_U.determinant())/2.0; // 1/2 is the area of the parametric space
        dA_list.push_back(dA);
        IU_list.push_back(I_U);
        ID_list.push_back(I_D);
        
        for(int k=0;k<3;k++)
        {
            if(TT(i,k)!=-1)
            {
                index[(k+2)%3] = F( TT(i,k), (TTi(i,k)+2)%3 );
                adjacent_points[(k+2)%3] = V.row(index[(k+2)%3]);
                adjacent_points_U[(k+2)%3] = V0.row(index[(k+2)%3]);
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
        // GeoFeature::calculate_second_fundamental_form(V0.row(F(i,0)), V0.row(F(i,1)), V0.row(F(i,2)), adjacent_points_U[0], adjacent_points_U[1], adjacent_points_U[2], real_pts, II_U);

        // IIU_list.push_back(II_U);
        IID_list.push_back(II_D);

    }
    
    // compute the derivative for the streching energy, which serves as the first term of the LHS equation
    for(int i=0;i<V.rows();i++)
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
//                A.coeffRef(3 * i + k,f) += 1.0 / 2.0 * thickness * dA_list[f] *
//                (0.5 * alpha * Q.trace() * dQ[k].trace() +
//                 beta * (Q * dQ[k]).trace());
                A(3 * i + k,f) += 1.0 / 2.0  * thickness * dA_list[f] *
                                (0.5 * alpha * Q.trace() * dQ[k].trace() +
                                 beta * (Q * dQ[k]).trace()) + 1.0 / 2.0 * thickness * dA_list[f] * (alpha+beta) * dQ[k].trace();
                
                b(3 * i + k) += 1.0 / 2.0 * thickness * dA_list[f] * (alpha+beta) * dQ[k].trace();
                
            }
        }
    }
    
    // compute the derivative for the bending energy, which serves as the second term of the LHS equation
    for(int i=0;i<V.rows();i++)
    {
        for(int j=0;j<VF[i].size();j++) // Consider the faces which contains the vertex
        {
            Eigen::Matrix2d Q;
            std::vector<Eigen::Matrix2d> dII,dQ;
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
            //Q = IU_list[f].inverse()*(IID_list[f]-IIU_list[f]);
            Q = IU_list[f].inverse()*IID_list[f];
            dQ.resize(3);
            for(int k=0;k<3;k++)
            {
                dQ[k] = IU_list[f].inverse()*dII[k];
//                A.coeffRef(3 * i + k,f) += 1.0 / 6.0 * pow(thickness,3.0) * dA_list[f] *
//                (0.5 * alpha * Q.trace() * dQ[k].trace() +
//                 beta * (Q * dQ[k]).trace());
                A(3 * i + k,f) += 1.0 / 6.0 * pow(thickness,3.0) * dA_list[f] *
                (0.5 * alpha * Q.trace() * dQ[k].trace() +
                 beta * (Q * dQ[k]).trace());

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
                //Q = IU_list[fa].inverse()*(IID_list[fa]-IIU_list[fa]);
                Q = IU_list[fa].inverse()*IID_list[fa];
                dQ.resize(3);
                for(int k=0;k<3;k++)
                {
                    dQ[k] = IU_list[fa].inverse()*dII[k];
//                    A.coeffRef(3 * i + k,fa) += 1.0 / 6.0 * pow(thickness,3.0) * dA_list[fa] *
//                    (0.5 * alpha * Q.trace() * dQ[k].trace() +
//                     beta * (Q * dQ[k]).trace());
                    A(3 * i + k,fa) += 1.0 / 6.0 * pow(thickness,3.0) * dA_list[fa] *
                    (0.5 * alpha * Q.trace() * dQ[k].trace() +
                     beta * (Q * dQ[k]).trace());

                }
            }
        }
    }
    //A.makeCompressed();
}

Eigen::VectorXd ComputeCoefficient::find_optimal_sol(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, double YonungsModulus, double PoissonRatio, double thickness, double M1, double M2)
{
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    construct_LS(V0, V, F, A, b, YonungsModulus, PoissonRatio, thickness);
    Eigen::MatrixXd G,CE,CI;
    Eigen::VectorXd g0,ce0,ci0,optimal_sol;
    G = A.transpose()*A;
    g0 = -b.transpose()*A;
    CE.resize(A.cols(), A.cols());
    CE.setZero();
    ce0.resize(A.cols());
    ce0.setZero();
    
    CI.resize(2*A.cols(), A.cols());
    CI.setZero();
    for(int j=0;j<CI.cols();j++)
    {
        CI(j,j)=1;
        CI(j+CI.cols(),j)=1;
    }
    ci0.resize(2*A.cols());
    for(int i=0;i<A.cols();i++)
    {
        ci0(i) = -M1;
        ci0(i+A.cols()) = M2;
    }
    optimal_sol.resize(A.cols());
    Eigen::solve_quadprog(G, g0, CE, ce0, CI, ci0, optimal_sol);
    std::cout<<"found omega!!"<<std::endl;
    return optimal_sol;
}

Eigen::VectorXd ComputeCoefficient::get_coefficient(Eigen::MatrixXd V0, Eigen::MatrixXd V, Eigen::MatrixXi F, double YonungsModulus, double PoissonRatio, double thickness, double M1, double M2)
{
    using namespace alglib;
    //Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    construct_LS(V0, V, F, A, b, YonungsModulus, PoissonRatio, thickness);
    
    // Use the optimization functions atated in alglib to solve this problem
    Eigen::MatrixXd Quad_A;
    Eigen::VectorXd rhs;
    
    Quad_A = A.transpose()*A;
    std::cout<<Quad_A.norm()<<std::endl;
    rhs = A.transpose()*b;

    
    //sparsematrix a;
    real_2d_array a;
    real_1d_array b1;
    real_1d_array x0;
    real_1d_array s;
    real_1d_array x;
    real_1d_array lb;
    real_1d_array ub;
    minqpstate state;
    minqpreport rep;

    int num_f = F.rows();
    b1.setlength(num_f);
    x0.setlength(num_f);
    lb.setlength(num_f);
    ub.setlength(num_f);
    s.setlength(num_f);

    for(int i=0;i<num_f;i++)
    {
        x0[i] = 1;
        lb[i] = M1;
        ub[i] = M2;
        s[i] = 1;
        b1[i] = rhs(i);
    }

    a.setlength(num_f, num_f);
    for(int i=0;i<num_f;i++)
    for(int j=0;j<num_f;j++)
    {
        a[i][j] = Quad_A(i,j);
    }
    // initialize sparsematrix structure
//    sparsecreate(num_f, num_f, 0, a);
//    for (int k=0; k<Quad_A.outerSize(); ++k)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(Quad_A,k); it; ++it)
//        {
////        it.value(); // value
////        it.row();   // row index
////        it.col();   // col index = k
//            sparseset(a, it.row(), it.col(), it.value());
//        }

    // create solver, set quadratic/linear terms, constraints
    minqpcreate(num_f, state);
    // minqpsetquadratictermsparse(state, a, true);
    minqpsetquadraticterm(state, a);
    minqpsetlinearterm(state, b1);
    minqpsetstartingpoint(state, x0);
    minqpsetbc(state, lb, ub);

    // Set scale of the parameters.
    // It is strongly recommended that you set scale of your variables.
    // Knowing their scales is essential for evaluation of stopping criteria
    // and for preconditioning of the algorithm steps.
    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
    //
    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
    //       which automatically determines variable scales.
    minqpsetscale(state, s);

    //
    // Solve problem with BLEIC-based QP solver.
    //
    // This solver is intended for problems with moderate (up to 50) number
    // of general linear constraints and unlimited number of box constraints.
    // It also supports sparse problems.
    //
    // Default stopping criteria are used.
    //
    // minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
    minqpoptimize(state);
    minqpresults(state, x, rep);
//    OUTPUT PARAMETERS:
//    X       -   array[0..N-1], solution.
//    This array is allocated and initialized only when
//    Rep.TerminationType parameter is positive (success).
//    Rep     -   optimization report. You should check Rep.TerminationType,
//    which contains completion code, and you may check  another
//    fields which contain another information  about  algorithm
//    functioning.
//
//    Failure codes returned by algorithm are:
//    * -9    failure of the automatic scale evaluation:  one of
//    the diagonal elements of  the  quadratic  term  is
//    non-positive.  Specify variable scales manually!
//    * -5    inappropriate solver was used:
//    * QuickQP solver for problem with  general  linear
//    constraints
//    * -4    BLEIC-QP/QuickQP   solver    found   unconstrained
//    direction  of   negative  curvature  (function  is
//                                          unbounded from below even under constraints),   no
//    meaningful minimum can be found.
//    * -3    inconsistent constraints (or maybe  feasible point
//                                      is too  hard  to  find).  If  you  are  sure  that
//    constraints are feasible, try to restart optimizer
//    with better initial approximation.
//
//    Completion codes specific for Cholesky algorithm:
//    *  4   successful completion
//
//    Completion codes specific for BLEIC/QuickQP algorithms:
//    *  1   relative function improvement is no more than EpsF.
//    *  2   scaled step is no more than EpsX.
//    *  4   scaled gradient norm is no more than EpsG.
//    *  5   MaxIts steps was taken
    Eigen::VectorXd coef(num_f);
    for(int i=0;i<num_f;i++)
    {
        coef(i) = x[i];
    }

    printf("%d\n", int(rep.terminationtype));
    return coef;
}

void ComputeCoefficient::test()
{
    Eigen::MatrixXd V0,V;
    Eigen::MatrixXi F0,F;
    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/DrapedRect/3876_triangles/draped_rect_geometry.obj", V0, F0);
    igl::readOBJ("/Users/chenzhen/UT/Research/Results/mine_simple_3.obj", V, F);
    
    Eigen::VectorXd omega_vec(F0.rows());
    
    // Calculate the energy
    for(int i=0;i<F0.rows();i++)
    {
        Eigen::Vector3d barycenter = (V0.row(F(i,0)) + V0.row(F(i,1)) + V0.row(F(i,2)))/3.0;
        double omega = 16.0*std::pow(0.025,4.0) / std::pow(4*0.025*0.025 + barycenter(0)*barycenter(0) + barycenter(1)*barycenter(1),2.0);
        omega_vec(i) = 1.0/omega;
    }
    double thickness = 1e-4;
    double Youngsmodulus = 1.0e4;
    double PoissonRatio = 0.3;
    Eigen::VectorXd coef = find_optimal_sol(V0,V,F,Youngsmodulus,PoissonRatio,thickness,1e-7,3);
    
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    construct_LS(V0, V, F, A, b, Youngsmodulus, PoissonRatio, thickness);
    //Eigen::VectorXd coef_1 = A.fullPivHouseholderQr().solve(b);
    std::cout<<coef<<std::endl;
    std::cout<<omega_vec<<std::endl;
    std::cout<<(A*coef-b).norm()<<std::endl;
    std::cout<<(A*omega_vec-b).norm()<<std::endl;
    std::cout<<(coef-omega_vec).norm()<<std::endl;
    std::cout<<"finished!"<<std::endl;
}
