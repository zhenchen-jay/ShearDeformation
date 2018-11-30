#include "LSQR.h"
#include "external/LSQR/lsqrDense.h"
#include "external/LSQR/lsmrDense.h"
#include <iostream>
#include <memory>
using namespace Eigen;
#ifndef MAX_ITR_TIME
#define MAX_ITR_TIME 1e1
#endif

VectorXd LSQR::LSQR_solver(SparseMatrix<double> A, VectorXd b, VectorXd x_ini, double eps)
{
    auto solver = std::make_unique<lsmrDense>();
    solver->SetOutputStream( std::cout );
    solver->SetEpsilon( eps );
    solver->SetDamp( 0.0 );
    solver->SetMaximumNumberOfIterations( 1000 );
    solver->SetToleranceA( 1e-16 );
    solver->SetToleranceB( 1e-16 );
    solver->SetUpperLimitOnConditional( 1.0 / ( 10 * sqrt( eps ) ) );
    //solver->SetStandardErrorEstimatesFlag( true );
    
    int rows = A.rows();
    int cols = A.cols();
    
   // double se[cols];
   // solver->SetStandardErrorEstimates( se );
    
    double bb[rows];
    double xx[cols];
    for(int i=0;i<rows;i++)
    {
        bb[i] = b(i);
    }
    
    double **AA;
    AA = new double *[rows];
    for(int i=0;i<rows;i++)
    {
        AA[i] = new double [cols];
    }
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
        {
            AA[i][j] = A.coeff(i, j);
        }
    solver->SetMatrix(AA);
    solver->Solve(rows, cols, bb, xx);
    VectorXd x(cols);
    x.setZero();
    for(int i=0;i<cols;i++)
    {
        x(i) = xx[i];
    }
    for(int i=0;i<rows;i++)
    {
        delete [] AA[i];
    }
    return x;
//    /******************************
//    * Initialize
//    ******************************/
//
//    VectorXd x = x_ini;
//    VectorXd x_old = x_ini;
//
//    double beta = (b - A * x_old).norm();
//    VectorXd u = (b - A * x_old) / beta;
//    std::cout<<u<<std::endl;
//    VectorXd ATu = A.transpose() * u;
//    std::cout<<ATu<<std::endl;
//    double alpha = ATu.norm();
//    VectorXd v = ATu/alpha;
//    VectorXd w = v;
//    double phi_bar = beta;
//    double rho_bar = alpha;
//
//    /***
//    * Variables for stopping criteria
//    ****/
//    double z = 0;
//    double cs2 = -1;
//    double sn2 = 0;
//    double ddnorm = 0;
//    double bnorm = beta;
//    double rnorm = beta;
//    double xnorm = 0;
//    double xxnorm = 0;
//    double Anorm = 0;
//    double Acond = 0;
//
//    int itr = 0;
//    while (itr < MAX_ITR_TIME){
//        /*************************************
//        * Continue the bidiagnolization
//        **************************************/
//        VectorXd rhs_beta = A *v - alpha * u;
//        std::cout<<rhs_beta<<std::endl;
//        beta = rhs_beta.norm();
//        u = rhs_beta / beta;
//
//        VectorXd rhs_alpha = A.transpose() * u  - beta * v;
//        alpha = rhs_alpha.norm();
//        v = rhs_alpha / alpha;
//
//        /*************************************
//        * Constract and apply next orthogonal transformation
//        **************************************/
//
//        double rho = sqrt(rho_bar * rho_bar + beta * beta);
//        double c = rho_bar / rho;
//        double s = beta / rho;
//        double theta = s * alpha;
//        rho_bar = -c* alpha;
//        double phi = c * phi_bar;
//        phi_bar = s*phi_bar;
//
//
//
//
//        /*************************************
//        * Test for convergence
//        **************************************/
//
//        double gambar = -cs2 *rho;
//        double rhs = phi - sn2 * rho * z;
//        double zbar = rhs / gambar;
//        xnorm = sqrt(xxnorm + zbar * zbar);
//        double gamma = sqrt(gambar* gambar + theta* theta);
//        cs2 = gambar / gamma;
//        sn2 = theta / gamma;
//        z = rhs / gamma;
//        xxnorm += z * z;
//
//
//        VectorXd rhow = (1 / rho) * w;
//        ddnorm = ddnorm + rhow.norm() * rhow.norm();
//        Anorm = sqrt(Anorm * Anorm + alpha * alpha + beta * beta);
//        Acond = Anorm + sqrt(ddnorm);
//        rnorm = phi_bar;
//        double Arnorm = alpha * abs(s * phi);
//        double test1 = rnorm / bnorm;
//        double test2 = 0;
//        double test3 = 0;
//        if (Anorm == 0 || rnorm == 0){
//            test2 = std::numeric_limits<double>::max();
//        }
//        else{
//            test2 = Arnorm / (Anorm * rnorm);
//        }
//        if (Acond == 0){
//            test3 = std::numeric_limits<double>::max();
//        }
//        else{
//            test3 = 1 / Acond;
//        }
//        double rtol = eps + eps * Anorm * xnorm / bnorm;
//
//
//        itr++;
//        if (test1 <= rtol || test2 <= eps || test3 <= eps){
//            if(test1 <= rtol)
//            {
//                std::cout<<"The relative error of r = b - Ax is small, itereate "<<itr<<" times"<<std::endl;
//                std::cout<<rnorm<<" "<<bnorm<<" "<<rnorm/bnorm<<std::endl;
//            }
//            else if(test2 <= eps)
//            {
//                std::cout<<"The relative error of Ar = A(b - Ax) is small, itereate "<<itr<<" times"<<std::endl;
//                std::cout<<(A.transpose()*(b-A*x)).norm()<<std::endl;
//            }
//            else
//            {
//                std::cout<<"The condotion number of the approximate matrix reach the bar, itereate "<<itr<<" times"<<std::endl;
//            }
//            break;
//        }
//
//        /*************************************
//        * Update x, w
//        **************************************/
//        x = x_old + (phi / rho) * w;
//        w = v - (theta / rho) * w;
//
//
//        // update history of x
//        x_old = x;
//
//    }
//    if(itr == MAX_ITR_TIME)
//        std::cout<<"terminate with the maximun iteration times"<<std::endl;
//    return x;
}

VectorXd LSQR::LSQR_solver_box_constraints(SparseMatrix<double> A, VectorXd b, VectorXd x_ini, VectorXd lx, VectorXd ux, double eps, double beta)
{
    VectorXd x = LSQR_solver(A, b,x_ini,eps);
    std::cout<<x<<std::endl<<std::endl;
    Projection(x, lx, ux);
    VectorXd r0 = A*x-b;
    VectorXd r = A.transpose()*(A*x-b);
    VectorXd lambda_l(x.size()), Al(x.size());
    VectorXd lambda_u(x.size()), Au(x.size());
    lambda_l.setZero();
    Al.setZero();
    lambda_u.setZero();
    Au.setZero();
    SparseMatrix<double> D(A.cols(),A.cols());
    D.setZero();
    
    VectorXd z_ini(x.size());
    z_ini.setZero();
    int itr = 0;
    
    while(r.norm() >= beta*eps && itr < MAX_ITR_TIME)
    {
        for(int i=0;i<x.size();i++)
        {
            if(x(i) == lx(i))
            {
                lambda_l(i) = r(i);
                if(r(i) > 0)
                    Al(i) = 1;
            }
            if(x(i) == ux(i))
            {
                lambda_u(i) = r(i);
                if(r(i) < 0)
                    Au(i) = 1;
            }
            D.coeffRef(i, i) = (1 - Al(i))*(1-Au(i));
        }
        
        VectorXd z = LSQR_solver(A*D, -r0, z_ini, beta*eps);
        x = x + D*z;
        Projection(x, lx, ux);
        r0 = A*x - b;
        r = A.transpose() * r0;
        itr++;
        std::cout<<"Iteration time: "<<itr<<" "<<z.norm()<<std::endl;
        std::cout<<x<<std::endl;
        std::cout<<"The norm of the gradient is: "<<r.norm()<<std::endl;
    }
    if(itr == MAX_ITR_TIME)
    {
        std::cout<<"Maximun iteration times reached!!"<<std::endl;
        std::cout<<"The norm of gradient at optimal solution is: "<<(A.transpose()*(A*x-b)).norm()<<std::endl;
    }
    else
    {
        std::cout<<"Terminate when the norm of gradient is small: "<<(A.transpose()*(A*x-b)).norm()<<std::endl;
        std::cout<<"Iteration times: "<<itr<<std::endl;
    }
    return x;
}

void LSQR::Projection(VectorXd &x, VectorXd lx, VectorXd ux)
{
    if(x.size() != lx.size() || x.size() != ux.size())
    {
        std::cout<<"The size of the box constraints don't match the size of the vector"<<std::endl;
        std::cout<<x.size()<<" "<<lx.size()<<" "<<ux.size()<<std::endl;
        return;
    }
    for(int i=0;i<x.size();i++)
    {
        if(x(i) < lx(i))
            x(i) = lx(i);
        else if(x(i) > ux(i))
            x(i) = ux(i);
    }
    
}
