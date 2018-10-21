//
// Created by Zhen Chen on 10/20/18.
//

#ifndef PROJECT_ELLIPTIC_FUNCTIONS_H
#define PROJECT_ELLIPTIC_FUNCTIONS_H
#include <complex.h>

// This function is based on the source code provided on "http://lasp.colorado.edu/cism/CISM_DX/code/CISM_DX-0.50/required_packages/octave-forge/main/specfun/ellipj.cc"

void jacobi_sn_cn_dn(double u, double m, double &sn, double &cn, double &dn, double &err);

void jacobi_sn_cn_dn(std::complex<double> u, double m, std::complex<double> &sn, std::complex<double> &cn, std::complex<double> &dn, double &err);

#endif //PROJECT_ELLIPTIC_FUNCTIONS_H
