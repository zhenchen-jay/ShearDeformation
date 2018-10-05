//
// Created by Zhen Chen on 10/3/18.
//

#include "ShellSimulation.h"
#include "ShellEnergy.h"
#include "GravityEnergy.h"

void callback_func(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df, void* ptr)
{
    ShellSimulation* callback = (ShellSimulation*) ptr;
    callback->energy_func_grad(x,f,df);
}


void ShellSimulation::compute_deformed_surface(Eigen::MatrixXd &V0, Eigen::MatrixXi &F0, double PossionRatio, double YoungsModulus, double thickness, double ratio)
{
    VU=V0;
    for(int i=0;i<V0.rows();i++)
    {
        if(abs(V0(i,0)-0.025)<1e-6 || abs(V0(i,0)+0.025)<1e-6)
        {
            Eigen::Vector3d point(ratio* V0(i,0),V0(i,1),V0(i,2));
            p_fixed.push_back(point);
            p_fixed_index.push_back(i);
            V0.row(i) = point.transpose();
        }
    }
    F=F0;
    _PossionRation = PossionRatio;
    _YoungsModulus = YoungsModulus;
    _thickness = thickness;
    _is_initialized = true;
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
    auto op_shell = std::make_unique<ShellEnergy>();
    double E = 0;
    Eigen::VectorXd dE;
    op_shell->streching_energy(V0, VU, F0, 1e5, 0.3, 1, E, dE);
}

void ShellSimulation::energy_func_grad(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df)
{
    if(x.length()%3!=0)
    {
        std::cout<<"Error!"<<std::endl;
        return;
    }
    if(!_is_initialized)
    {
        std::cout<<"Please do initialization first!"<<std::endl;
        return;
    }
    //convert x to vertice format
    int v_num = int(x.length()/3);
    Eigen::MatrixXd V(v_num,3);
    for(int i=0;i<v_num;i++)
    {
        V.row(i)<<x[3*i],x[3*i+1],x[3*i+2];
    }
    auto op_shell = std::make_unique<ShellEnergy>();
    auto op_gravity = std::make_unique<GravityEnergy>();
    double E_shell(0), E_gravity(0);
    Eigen::VectorXd diff_f_shell,diff_f_gravity;
    op_shell->streching_energy(V,VU,F,_YoungsModulus,_PossionRation,_thickness,E_shell,diff_f_shell);
    op_gravity->gravity_energy(V, 1, E_gravity, diff_f_gravity);
    f = E_shell;
    for(int i=0;i<df.length();i++)
    {
        df[i] = diff_f_shell(i);
    }
    for(int i=0;i<p_fixed_index.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            df[3*p_fixed_index[i]+j]=0;
        }
    }
}


