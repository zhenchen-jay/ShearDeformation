// #include <ifopt/variable_set.h>
// #include <ifopt/constraint_set.h>
// #include <ifopt/cost_term.h>
// #include <iostream>

// #ifndef EPS_BOUND
// #define EPS_BOUND 1e-10
// #endif

// namespace ifopt
// {
//     class optVariables: public VariableSet
//     {
//     public:
//         optVariables(int num, const std::string& name) : VariableSet(num, name)
//         {
//             _x.resize(num);
//             _x.setZero();
//             if(num%3 != 0)
//             {
//                 std::cout<<"Wrong variable size"<<std::endl;
//             }
//             else
//             {
//                 int num_f =  num/3;
//                 for(int i=0;i<num_f;i++)
//                 {
//                     _x(3*i) = 1;
//                     _x(3*i+2) = 1;
//                 }
//             }
//         }
        
//         void SetVariables(const VectorXd& x) override
//         {
//             _x = x;
//         }
        
//         VectorXd GetValues() const override
//         {
//             return _x;
//         }
        
//         VecBound GetBounds() const override
//         {
//             VecBound bounds(GetRows());
//             int num_f = GetRows()/3;
//             for(int i=0;i<num_f;i++)
//             {
//                 bounds.at(3*i) = Bounds(EPS_BOUND,inf);
//                 bounds.at(3*i+1) = NoBound;
//                 bounds.at(3*i+2) = NoBound;
//             }
//             return bounds;
//         }
        
//     private:
//         VectorXd _x;
//     };
    
//     class optConstraint: public ConstraintSet
//     {
//     public:
//         optConstraint(int num, const std::string& name) : ConstraintSet(num, name){}
        
//         VectorXd GetValues() const override
//         {
//             VectorXd g(GetRows());
//             VectorXd x = GetVariables()->GetComponent("var_set")->GetValues();
//             if(g.size()*3!=x.size())
//             {
//                 std::cout<<"Invalid constraints number"<<std::endl;
//             }
//             else
//             {
//                 for(int i=0;i<g.size();i++)
//                 {
//                     g(i) = x(3*i)*x(3*i+2) - x(3*i+1)*x(3*i+1);
//                 }
//             }
//             return g;
//         }
        
//         VecBound GetBounds() const override
//         {
//             VecBound b(GetRows());
//             for(int i=0;i<GetRows();i++)
//             {
//                 b.at(i) = Bounds(EPS_BOUND,inf);
//             }
//             return b;
//         }
        
//         void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override
//         {
//             if(var_set == "var_set")
//             {
//                 VectorXd x = GetVariables()->GetComponent("var_set")->GetValues();
// //                int num_f = GetRows();
//                 for(int i=0;i<GetRows();i++)
//                 {
//                     jac_block.coeffRef(0, 3*i) = x(3*i+2);
//                     jac_block.coeffRef(0, 3*i+1) = -2*x(3*i+1);
//                     jac_block.coeffRef(0, 3*i+2) = x(3*i);
//                 }
//             }
//         }
//     };
    
//     class optCost: public CostTerm
//     {
//     public:
//         optCost(Eigen::MatrixXd V1, Eigen::MatrixXi F1, double YoungsModulus, double PoissonRatio, double thickness) : optCost("cost_term1")
//         {
//             _is_initialized = false;
//             set_up(V1, F1, YoungsModulus, PoissonRatio, thickness);
//             _is_initialized = true;
//         }
//         optCost(const std::string& name) : CostTerm(name){}
        
//         double GetCost() const override;
//         void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
        
//     private:
//         void compute_derivative_inv_mat(Eigen::Matrix2d A, std::vector<Eigen::Matrix2d > &dA) const;
//         void compute_derivative_sqrt_det(Eigen::Matrix2d A, std::vector<double> &diff_sqrt_det) const;
//         void set_up(Eigen::MatrixXd VD, Eigen::MatrixXi F0, double YonungsModulus, double PoissonRatio, double thickness);
        
        
//     private:
//         bool _is_initialized;
//         std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
//         std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
//         Eigen::MatrixXi TT;             // The matrix which stores the information of adjacent faces
//         Eigen::MatrixXi TTi;
//         Eigen::MatrixXd V;
//         Eigen::MatrixXi F;
        
//         std::vector<Eigen::Matrix2d> ID_list;
//         std::vector<Eigen::Matrix2d> IID_list;
//         //
//         std::vector<std::vector<Eigen::Matrix2d > > dID_list; // The list to store the derivative of the first fundamental form. The order is determined by VF
//         std::vector<std::vector<Eigen::Matrix2d > > dIID_list_neighbor;// The list to store the derivative of the first fundamental form. The order is determined by VF
//         std::vector<std::vector<Eigen::Matrix2d > > dIID_list_2_neighbor; // The list to store the derivative of the first fundamental form. The order is determined by TT
        
//         double _PoissonsRatio;
//         double _YoungsModulus;
//         double _thickness;
//         double _alpha;
//         double _beta;
        
//     };
    
// }
