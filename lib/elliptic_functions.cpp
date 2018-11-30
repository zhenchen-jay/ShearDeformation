//
// Created by Zhen Chen on 10/20/18.
//
#if 0
#include "elliptic_functions.h"
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <iostream>
#include <cmath>

#ifndef M_EPS
#define M_EPS 2.220446049e-16
#endif

#ifndef NMAX
#define NMAX 16
#endif

void jacobi_sn_cn_dn(double u, double m, double &sn, double &cn, double &dn, double &err)
{
    double sqrt_eps, m1, t=0., si_u, co_u, se_u, ta_u, b, c[NMAX], a[NMAX], phi;
    int n, Nn, ii;
    
    if (m < 0. || m > 1.)
    {
        std::cout<<"Invalid input for m, which should be between 0 and 1"<<std::endl;
        return;
    }
    sqrt_eps = sqrt(M_EPS);
    if (m < sqrt_eps)
    {
        /*  # For small m, ( Abramowitz and Stegun, Section 16.13 ) */
        /*{{{*/
        si_u = sin(u);
        co_u = cos(u);
        t = 0.25*m*(u-si_u*co_u);
        sn = si_u - t * co_u;
        cn = co_u + t * si_u;
        dn = 1.0 - 0.5*m*si_u*si_u;
        /*}}}*/
    }
    else if ( (1.0 - m) < sqrt_eps )
    {
        /*  For m1 = (1-m) small ( Abramowitz and Stegun, Section 16.15 ) */
        /*{{{*/
        m1 = 1.0-m;
        si_u = sinh(u);
        co_u = cosh(u);
        ta_u = tanh(u);
        se_u = 1.0/co_u;
        sn = ta_u + 0.25*m1*(si_u*co_u-u)*se_u*se_u;
        cn = se_u - 0.25*m1*(si_u*co_u-u)*ta_u*se_u;
        dn = se_u + 0.25*m1*(si_u*co_u+u)*ta_u*se_u;
        /*}}}*/
    }
    else
    {
        a[0] = 1.0;
        b    = sqrt(1.0-m);
        c[0] = sqrt(m);
        for (n = 1; n<NMAX; ++n)
        {
            a[n] = (a[n-1]+b)/2;
            c[n] = (a[n-1]-b)/2;
            b = sqrt(a[n-1]*b);
            if ( c[n]/a[n] < M_EPS) break;
        }
        if ( n >= NMAX-1)
        {
            std::cout<<"Please increase your iteration times"<<std::endl;
            err = 1.;
            return;
        }
        Nn = n;
        for ( ii = 1;  n>0;    ii = ii*2, --n) ; // pow(2, Nn)
        phi = ii*a[Nn]*u;
        for ( n = Nn; n > 0; --n)
        {
            t = phi;
            phi = (asin((c[n]/a[n])* sin(phi))+phi)/2.;
        }
        sn = sin(phi);
        cn = cos(phi);
        dn = cn/cos(t-phi);
    }
    return;
}

void jacobi_sn_cn_dn(std::complex<double> u, double m, std::complex<double> &sn, std::complex<double> &cn, std::complex<double> &dn, double &err)
{
    using namespace std::complex_literals;
    double m1 = 1.-m, ss1, cc1, dd1;
    jacobi_sn_cn_dn(u.imag(), m1, ss1, cc1, dd1, err);
    if ( abs(u.real()) < 1e-10)
    {
        /* u is pure imag: Jacoby imag. transf. */
        /*{{{*/
        sn = ss1/cc1*1i;
        cn = 1/cc1;         //    cn.imag = 0.;
        dn = dd1/cc1;       //    dn.imag = 0.;
    }
    else
    {
        /* u is generic complex */
        double ss, cc, dd, ddd;
        
        jacobi_sn_cn_dn( u.real(), m, ss, cc, dd, err);
        ddd = cc1*cc1 + m*ss*ss*ss1*ss1;
        sn = ss*dd1/ddd + cc*dd*ss1*cc1/ddd*1i;
        cn = cc*cc1/ddd - ss*dd*ss1*dd1/ddd*1i;
        dn = dd*cc1*dd1/ddd - m*ss*cc*ss1/ddd*1i;
    }
    return;
}
#endif
