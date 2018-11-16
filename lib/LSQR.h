#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace LSQR
{
    using namespace Eigen;
    VectorXd LSQR_solver(SparseMatrix<double> A, VectorXd b, VectorXd x_ini, double eps);
    VectorXd LSQR_solver_box_constraints(SparseMatrix<double> A, VectorXd b, VectorXd x_ini, VectorXd lx, VectorXd ux, double eps, double beta);
    void Projection(VectorXd &x, VectorXd lx, VectorXd ux);
}

