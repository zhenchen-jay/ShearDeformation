#include <igl/vertex_triangle_adjacency.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <iomanip>
#include "simulation.h"

#ifndef M_FORCE
#define M_FORCE 10
#endif

int times = 0;

void callback_func(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df, void* ptr)
{
    Simulation* callback = (Simulation*) ptr;
    callback->energy_func_grad(x,f,df);
}

void Simulation::test(Eigen::MatrixXd &V, Eigen::MatrixXi F0,double nu, double E, double t, double ratio)
{
    // Initialize
    //V_U = ratio*V_U_record;
    if(!_is_initialized)
        initialize(V, F0, nu, E, t);
    else
    {
        V_U=V;
    }
    
    // get the IU_list
    IU_list.resize(F.rows());
    dA_list.resize(F.rows());
    for(int i=0;i<F.rows();i++)
    {
        calculate_fundamental_form1(V_U.row(F(i,0)),V_U.row(F(i,1)),V_U.row(F(i,2)),IU_list[i]);
        dA_list[i] = sqrt(IU_list[i](0,0)*IU_list[i](1,1)-IU_list[i](0,1)*IU_list[i](1,0));
    }
    
    p_fixed_index.clear();
    p_fixed.clear();
    
    for(int i=0;i<V_U_record.rows();i++)
    {
        if(abs(V_U_record(i,0)-0.025)<1e-6 || abs(V_U_record(i,0)+0.025)<1e-6)
        {
            Eigen::Vector3d point(ratio* V_U_record(i,0),V_U_record(i,1),V_U_record(i,2));
            p_fixed.push_back(point);
            p_fixed_index.push_back(i);
        }
    }
    
    for(int i=0;i<p_fixed_index.size();i++)
    {
        V.row(p_fixed_index[i])=p_fixed[i];
    }
    //std::cout<<membrane_energy_default(V)<<std::endl;
    //std::cout<<diff_membrane_energy_default(V).norm()<<std::endl;
    
    // Simulation
    alglib::real_1d_array x;
    double epsg = 1e-15;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0.1;
    alglib::ae_int_t maxits = 30000;
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    x.setlength(3*V.rows());
    for(int i=0;i<V.rows();i++)
        for(int j=0;j<3;j++)
        {
            x[3*i+j] = V(i,j);
        }

    // Solve the optimization problem using alglib
    alglib::minlbfgscreate(1, x, state);
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    alglib::minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, callback_func, NULL, this);
    alglib::minlbfgsresults(state, x, rep);
    for(int i=0;i<V.rows();i++)
    {
        V.row(i)<<x[3*i],x[3*i+1],x[3*i+2];
    }
}

void Simulation::test_function(Eigen::VectorXd v, Eigen::MatrixXd V, double eps)
{
    V(0,1) = V(0,1) + 0.1*V(0,1);
    V(1,1) = V(1,1) + 0.1*V(1,1);
    
    //Eigen::VectorXd diff_external = diff_external_energy(V);
    Eigen::VectorXd diff_membrane = diff_membrane_energy_default(V);
    double grad_results = (diff_membrane).dot(v);
    double membrane_E1 = membrane_energy_default(V);
    double energe_result1 = membrane_E1;
    int v_num = V.rows();
    Eigen::VectorXd x(3*v_num);
    if (v.size()!=x.size())
    {
        std::cout<<"Wrong input"<<std::endl;
        return;
    }
    for (int i=0;i<v_num; i++)
    {
        x.segment(3*i, 3) = V.row(i).transpose();
    }
    x = x+eps*v;
    for (int i=0; i<v_num; i++)
    {
        V.row(i) = x.segment(3*i, 3).transpose();
    }
    double membrane_E2 = membrane_energy_default(V);
    double energe_result2 = membrane_E2;
    
    
    std::cout<<"The power to the epsilon: "<<int(log10(eps))<<std::endl;
    std::cout<<std::fixed<<std::setprecision(0)<<v.transpose()<<std::endl;
    std::cout.setf(std::ios::fixed);
    std::cout<<std::setprecision(20);
    std::cout<<"membrane energy: "<<"Original: "<<membrane_E1<<" "<<"After disturbance: "<<membrane_E2<<std::endl;
    std::cout<<"gradient: "<<grad_results<<std::endl;
    double diff_results = (energe_result2-energe_result1)/eps;
    std::cout<<"Difference of Energy values: "<<diff_results<<std::endl;
    std::cout<<"Error: "<<abs(diff_results-grad_results)<<std::endl<<std::endl;
}

void Simulation::compute_deformed_surface(Eigen::MatrixXd &V0, Eigen::MatrixXi &F0, double nu, double E, double t)
{
    initialize(V0, F0, nu, E, t);

    // Simulation
    alglib::real_1d_array x;
    double epsg = 0;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0.1;
    alglib::ae_int_t maxits = 30000;
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    x.setlength(3*V0.rows());
    for(int i=0;i<V0.rows();i++)
        for(int j=0;j<3;j++)
        {
            x[3*i+j] = V0(i,j);
        }
    
    // Solve the optimization problem using alglib
    alglib::minlbfgscreate(1, x, state);
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    alglib::minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, callback_func, NULL, this);
    alglib::minlbfgsresults(state, x, rep);
    for(int i=0;i<V0.rows();i++)
    {
        V0.row(i)<<x[3*i],x[3*i+1],x[3*i+2];
    }
}

void Simulation::initialize(Eigen::MatrixXd V0, Eigen::MatrixXi F0, double nu, double E, double t)
{
    if(_is_initialized)
    {
        std::cout<<"The data have been initialized!"<<std::endl;
        return;
    }
    // Initialize
    V_U = V0;
    V_U_record = V0;
    F = F0;
    _A = 1.0/2;
    _alpha = E*nu/((1.0+nu)*(1.0-2.0*nu));
    _beta = E/(2.0*(1+nu));
    //_M = E/(1.0-nu*nu);
    _t = t;
    _nu = nu;
    calculate_coeff_mat(e0,e1,e2,_q,_nu);
    
    // set fixed points
    p_fixed.push_back(V0.row(2));
    p_fixed_index.push_back(2);
    p_fixed.push_back(V0.row(3));
    p_fixed_index.push_back(3);
    
    // Parametric plane
    e0<<1,0;
    e1<<-1,1;
    e2<<0,-1;
    
    // get the information of vertex-triangle adjacency
    igl::vertex_triangle_adjacency(V0.rows(), F0, VF, VFi);
    
    // get the IU_list
    IU_list.resize(F.rows());
    dA_list.resize(F.rows());
    for(int i=0;i<F.rows();i++)
    {
        calculate_fundamental_form1(V0.row(F(i,0)),V0.row(F(i,1)),V0.row(F(i,2)),IU_list[i]);
        dA_list[i] = sqrt(IU_list[i](0,0)*IU_list[i](1,1)-IU_list[i](0,1)*IU_list[i](1,0));
    }
    _is_initialized = true;
}

void Simulation::energy_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df)
{
    if(x.length()%3!=0)
    {
        std::cout<<"Error!"<<std::endl;
    }
    //convert x to vertice format
    int v_num = int(x.length()/3);
    Eigen::MatrixXd V(v_num,3);
    for(int i=0;i<v_num;i++)
    {
        V.row(i)<<x[3*i],x[3*i+1],x[3*i+2];
    }
    //f = membrane_energy_default(V) + external_energy(V);
    f = membrane_energy_default(V);
    //Eigen::VectorXd diff_f = diff_membrane_energy_default(V) + diff_external_energy(V);
    Eigen::VectorXd diff_f = diff_membrane_energy_default(V);
    for(int i=0;i<diff_f.size();i++)
    {
        df[i] = diff_f(i);
    }
}

//double Simulation::membrane_energy(Eigen::MatrixXd V)
//{
//    if(_is_initialized!=true)
//    {
//        std::cout<<"Please call initialize function first"<<std::endl;
//        return 0;
//    }
//    double E = 0;
//    for(int index_i = 0;index_i<F.rows();index_i++)
//    {
//        /*  Calculate the Matrix for the first fundamemtal form
//         The parametric plane is set as span{(0,0),(1,0),(0,1)}*/
//        Eigen::Matrix2d I_D;
//        calculate_fundamental_form1(V.row(F(index_i,0)),V.row(F(index_i,1)),V.row(F(index_i,2)),I_D);
//        Eigen::Vector3d Q;
//        Q << calculate_QI(IU_list[index_i],I_D,e0), calculate_QI(IU_list[index_i],I_D,e1), calculate_QI(IU_list[index_i],I_D,e2);
//        double test = 0;
//        for(int o=0;o<3;o++)
//            for(int i=0;i<3;i++)
//            {
//                E = E + 1.0 * _t * _M * dA_list[index_i] * Q(i) * Q(o) * _q(i,o) / (256 * pow(_A,3.0));
//                test+= 1.0 * _t * _M * dA_list[index_i] * Q(i) * Q(o) * _q(i,o) / (256 * pow(_A,3.0));
//            }
//    }
//    return E;
//}

double Simulation::membrane_energy_default(Eigen::MatrixXd V)
{
    if(_is_initialized!=true)
    {
        std::cout<<"Please call initialize function first"<<std::endl;
        return 0;
    }
    double E = 0;
    for(int i=0;i<F.rows();i++)
    {
        Eigen::Matrix2d I_D;
        calculate_fundamental_form1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),I_D);
        Eigen::MatrixXd M = IU_list[i].inverse()*(I_D-IU_list[i]);
        E+=1.0 * _t *dA_list[i]*(_alpha*0.5*(M.trace()*M.trace())+_beta*(M*M).trace())/8;
        
    }
    std::cout<<E<<std::endl;
    return E;
}

//Eigen::VectorXd Simulation::diff_membrane_energy(Eigen::MatrixXd V)
//{
//    if(_is_initialized!=true)
//    {
//        std::cout<<"Please call initialize function first"<<std::endl;
//        return Eigen::VectorXd::Zero(1);
//    }
//    int v_num = V.rows();  // The number of the vertice
//    Eigen::VectorXd dE = Eigen::VectorXd::Zero(3*v_num);
//    // Compute dE
//    for (int i=0;i<v_num;i++)
//    {
//
//        for(int j=0;j<VF[i].size();j++) // Consider the faces which contains the vertex
//        {
//            Eigen::Matrix2d I_D, dI_x, dI_y, dI_z;
//            int f = VF[i][j];
//            int fi = VFi[i][j];
//            calculate_fundamental_form1(V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), I_D);
//            Eigen::Vector3d Q;
//            Eigen::Matrix3d dQ;
//            Q << calculate_QI(IU_list[f],I_D,e0), calculate_QI(IU_list[f],I_D,e1), calculate_QI(IU_list[f],I_D,e2);
//            dQ.row(0) = diff_QI(IU_list[f],V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), fi,e0);
//            dQ.row(1) = diff_QI(IU_list[f],V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), fi,e1);
//            dQ.row(2) = diff_QI(IU_list[f],V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), fi,e2);
//            for(int o=0;o<3;o++)
//                for(int k=0;k<3;k++)
//                {
//                    dE.segment(3*i,3) = dE.segment(3*i, 3) + 1.0 * _t * _M * dA_list[f] * _q(k,o) * (Q(o)*dQ.row(k)+Q(k)*dQ.row(o)).transpose() / (256 * pow(_A,3.0));
//                }
//        }
//    }
//
//    for(int i=0;i<p_fixed_index.size();i++)
//    {
//        dE(3*p_fixed_index[i])=0;
//        dE(3*p_fixed_index[i]+1)=0;
//        dE(3*p_fixed_index[i]+2)=0;
//    }
//    //std::cout<<dE.norm()<<std::endl;
//    return dE;
//}

Eigen::VectorXd Simulation::diff_membrane_energy_default(Eigen::MatrixXd V)
{
    if(_is_initialized!=true)
    {
        std::cout<<"Please call initialize function first"<<std::endl;
        return Eigen::VectorXd::Zero(1);
    }
    int v_num = V.rows();  // The number of the vertice
    Eigen::VectorXd dE = Eigen::VectorXd::Zero(3*v_num);
    // Compute dE
    for (int i=0;i<v_num;i++)
    {
        
        for(int j=0;j<VF[i].size();j++) // Consider the faces which contains the vertex
        {
            Eigen::Matrix2d I_D, dI_x, dI_y, dI_z,Q,dQ_x,dQ_y,dQ_z;
            int f = VF[i][j];
            int fi = VFi[i][j];
            calculate_fundamental_form1(V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), I_D);
            diff_fundamental_form1(V.row(F(f,0)), V.row(F(f,1)), V.row(F(f,2)), fi, dI_x, dI_y, dI_z);
            Q=IU_list[f].inverse()*(I_D-IU_list[f]);
            dQ_x =IU_list[f].inverse()*dI_x;
            dQ_y = IU_list[f].inverse()*dI_y;
            dQ_z = IU_list[f].inverse()*dI_z;
            dE(3*i) = dE(3*i) + 1.0*_t/4.0*dA_list[f]*(0.5*_alpha*Q.trace()*dQ_x.trace()+_beta*(Q*dQ_x).trace());
            dE(3*i+1) = dE(3*i+1) + 1.0*_t/4.0*dA_list[f]*(0.5*_alpha*Q.trace()*dQ_y.trace()+_beta*(Q*dQ_y).trace());
            dE(3*i+2) = dE(3*i+2) + 1.0*_t/4.0*dA_list[f]*(0.5*_alpha*Q.trace()*dQ_z.trace()+_beta*(Q*dQ_z).trace());
            
        }
    }
    
    for(int i=0;i<p_fixed_index.size();i++)
    {
        dE(3*p_fixed_index[i])=0;
        dE(3*p_fixed_index[i]+1)=0;
        dE(3*p_fixed_index[i]+2)=0;
    }
    std::cout<<dE.norm()<<std::endl;
    return dE;

}

double Simulation::external_energy(Eigen::MatrixXd V)
{
    double E = 0;
    for(int i=0;i<V.rows();i++)
    {
        if(V(i,1)+0.5>0)
            E = E + M_FORCE*(V(i,1)+0.5);
    }
//    for(int i=0;i<F.rows();i++)
//    {
//        if((V(F(i,0),1)+V(F(i,1),1)+V(F(i,2),1))/3.0+0.5>0)
//            E = E + 3.0*((V(F(i,0),1)+V(F(i,1),1)+V(F(i,2),1))/3.0+0.5)*M_FORCE;
//    }
    return E;
}

Eigen::VectorXd Simulation::diff_external_energy(Eigen::MatrixXd V)
{
    Eigen::VectorXd dE = Eigen::VectorXd::Zero(3*V.rows());
    for(int i=0;i<V.rows();i++)
    {
        if(V(i,1)+0.5>0)
            dE(3*i+1) = M_FORCE;
//        for(int j=0;j<VF[i].size();j++)
//        {
//            int f = VF[i][j];
//            if((V(F(f,0),1)+V(F(f,1),1)+V(F(f,2),1))/3.0+0.5>0)
//                dE(3*i+1) += M_FORCE;
//        }
    }
    for(int i=0;i<p_fixed_index.size();i++)
    {
        dE(3*p_fixed_index[i])=0;
        dE(3*p_fixed_index[i]+1)=0;
        dE(3*p_fixed_index[i]+2)=0;
    }
    return dE;
}

void Simulation::calculate_coeff_mat(Eigen::Vector2d v0, Eigen::Vector2d v1, Eigen::Vector2d v2, Eigen::Matrix3d &q, double nu)
{
    Eigen::MatrixXd per_mat(3,2);
    per_mat <<  0,-1,
                1,1,
                -1,0;
    for(int o=0;o<3;o++)
        for(int i=0;i<3;i++)
        {
            q(i,o) = get_coeff_mat_components(per_mat.row((i+1)%3), per_mat.row((i+2)%3), per_mat.row((o+1)%3), per_mat.row((o+2)%3), nu);
        }
}

double Simulation::get_coeff_mat_components(Eigen::Vector2d vj, Eigen::Vector2d vk, Eigen::Vector2d vp, Eigen::Vector2d vq, double nu)
{
    double p1 = 4*nu*vj.dot(vk)*vp.dot(vq);
    double p2 = 4*(1-nu)*(vj(0)*vk(0)*vp(0)*vq(0)+vj(1)*vk(1)*vp(1)*vq(1));
    double p3 = 2*(1-nu)*(vj(0)*vk(1)+vj(1)*vk(0))*(vp(0)*vq(1)+vp(1)*vq(0));
    return p1+p2+p3;
}

void Simulation::calculate_fundamental_form1(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, Eigen::Matrix2d &I)
{
    I(0,0) = (V1-V0).dot(V1-V0);
    I(1,0) = (V1-V0).dot(V2-V0);
    I(0,1) = (V1-V0).dot(V2-V0);
    I(1,1) = (V2-V0).dot(V2-V0);
}

void Simulation::diff_fundamental_form1(Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int flag, Eigen::Matrix2d &dI_x, Eigen::Matrix2d &dI_y, Eigen::Matrix2d &dI_z)
{
    if(flag!=0 && flag!=1 && flag!=2)
    {
        std::cout<<"The differential vertex should lie on the triangle"<<std::endl;
        return;
    }
    Eigen::Vector3d A11, A12, A21, A22;
    if(flag==0)
    {
        A11 = -2*(V1-V0);
        A12 = -(V1+V2-2*V0);
        A21 = -(V1+V2-2*V0);
        A22 = -2*(V2-V0);
    }
    else if(flag==1)
    {
        A11 = 2*(V1-V0);
        A12 = V2-V0;
        A21 = V2-V0;
        A22 << 0,0,0;
    }
    else
    {
        A11 << 0,0,0;
        A12 = V1-V0;
        A21 = V1-V0;
        A22 = 2*(V2-V0);
    }
    dI_x << A11(0),A12(0),
            A21(0),A22(0);
    dI_y << A11(1),A12(1),
            A21(1),A22(1);
    dI_z << A11(2),A12(2),
            A21(2),A22(2);

}

double Simulation::calculate_QI(Eigen::Matrix2d I_U, Eigen::Matrix2d I_D, Eigen::Vector2d v)
{
    return v.transpose()*I_U.inverse()*I_D*v-v.dot(v);
}

Eigen::Vector3d Simulation::diff_QI(Eigen::Matrix2d I_U, Eigen::Vector3d V0, Eigen::Vector3d V1, Eigen::Vector3d V2, int flag, Eigen::Vector2d e)
{
    Eigen::Matrix2d dI_x, dI_y, dI_z;
    diff_fundamental_form1(V0, V1, V2, flag, dI_x, dI_y, dI_z);
    Eigen::Vector3d dQ;
    dQ(0) = e.transpose()*I_U.inverse()*dI_x*e;
    dQ(1) = e.transpose()*I_U.inverse()*dI_y*e;
    dQ(2) = e.transpose()*I_U.inverse()*dI_z*e;
    return dQ;
}
