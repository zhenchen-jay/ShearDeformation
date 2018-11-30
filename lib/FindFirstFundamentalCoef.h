#ifndef FindFirstFundamentalCoef_h
#define FindFirstFundamentalCoef_h
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "external/alglib/stdafx.h"
#include "external/alglib/optimization.h"

class FindFirstFundamentalCoef
{
public:
    FindFirstFundamentalCoef(): _is_initialized(false),ID_list(0,Eigen::Matrix2d()),
    IID_list(0,Eigen::Matrix2d()), dID_index(0,std::vector<int>(0,0)),dIID_index(0,std::vector<int>(0,0)){}
    ~FindFirstFundamentalCoef(){}

public:
    void get_func_grad(Eigen::VectorXd &x, double &f, Eigen::VectorXd &df);
    void set_up(Eigen::MatrixXd VD, Eigen::MatrixXi F0, double YonungsModulus, double PoissonRatio, double thickness);
    
    void test_func_grad();

private:
    bool _is_initialized;
    std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
    std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
    Eigen::MatrixXi TT;             // The matrix which stores the information of adjacent faces
    Eigen::MatrixXi TTi;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    std::vector<Eigen::Matrix2d> ID_list;
    std::vector<Eigen::Matrix2d> IID_list;
    
    Eigen::SparseMatrix<double> dID_list;
    Eigen::SparseMatrix<double> dIID_list;
    
    std::vector<std::vector<int> > dID_index;
    std::vector<std::vector<int> > dIID_index;
    
    double _PoissonsRatio;
    double _YoungsModulus;
    double _thickness;
    double _alpha;
    double _beta;

};

#endif
