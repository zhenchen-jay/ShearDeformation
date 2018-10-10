////
//// Created by Zhen Chen on 9/17/18.
////
//
//#include <igl/readOFF.h>
//#include <igl/readOBJ.h>
//#include <igl/writeOBJ.h>
//#include <igl/opengl/glfw/Viewer.h>
//#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
//#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
//#include <igl/boundary_loop.h>
//#include <imgui/imgui.h>
////#include "lib/simulation.h"
//#include "lib/ShellEnergy.h"
//#include "lib/ShellSimulation.h"
//#include "lib/GeometryFeature.h"
//#include "lib/MeshConnection.h"
//#include <iostream>
//#include <string>
//#include <fstream>
//
//void write_boundary(Eigen::MatrixXd V, Eigen::MatrixXi F)
//{
//    std::vector<std::vector<int>> L;
//    igl::boundary_loop(F, L);
//    std::ofstream in;
//    in.open("boundary.txt",std::ios::trunc);
//    for(int i=0;i<L.size();i++)
//    {
//        in<<i<<std::endl; // The i-th boundary
//        for(int j=0;j<L[i].size();j++)
//        {
//            in<<V(L[i][j],0)<<" "<<V(L[i][j],1)<<" "<<V(L[i][j],2)<<std::endl;
//        }
//    }
//    in.close();
//}
//
//void test_grad()
//{
//    //auto op = std::make_unique<Simulation>();
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/plane.obj", V, F);
//    //op->initialize(V, F, 0.3, 1e5, 1);
//    for(int i=0;i<3*V.rows();i++)
//    {
//        if(i%3==2 || i == 7 || i==10 || i == 6 || i == 9)
//            continue;
//        std::cout<<std::endl;
//        Eigen::VectorXd v = Eigen::VectorXd::Zero(3*V.rows());
//        v(i)=1;
//        //for(int k = 3;k<10;k++)
//        //    op->test_function(v, V, pow(10,-k));
//    }
//}
//
//void test_sphere()
//{
//    Eigen::MatrixXd Vs;
//    Eigen::MatrixXi Fs;
//    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/sphere.obj", Vs, Fs);
//    std::vector<std::vector<int> > VF;  // list of lists of incident faces (adjacency list)
//    std::vector<std::vector<int> > VFi; // list of lists of index of incidence within incident faces listed in VF
//    Eigen::MatrixXi TT;             // The matrix which stores the information of adjacent faces
//    Eigen::MatrixXi TTi;
//    std::vector<Eigen::Matrix2d> I_list; // The list of the first fundamental form of the original surface
//    std::vector<Eigen::Matrix2d> II_list; // The list of the first fundamental form of the current surface
//    std::vector<int> index(3);
//    std::vector<Eigen::Vector3d> adjacent_points(3);
//    std::vector<bool> real_pts(3);
//
//    MeshConnection::vertex_triangle_adjacency(Vs.rows(), Fs, VF, VFi);
//    MeshConnection::triangle_triangle_adjacency(Fs, TT, TTi);
//
//    // Calculate the energy
//    for(int i=0;i<Fs.rows();i++)
//    {
//        Eigen::Matrix2d I, II;
//        GeoFeature::calculate_first_fundamental_form(Vs.row(Fs(i,0)),Vs.row(Fs(i,1)),Vs.row(Fs(i,2)),I);
//        I_list.push_back(I);
//        for(int k=0;k<3;k++)
//        {
//            if(TT(i,k)!=-1)
//            {
//                index[(k+2)%3] = Fs( TT(i,k), (TTi(i,k)+2)%3 );
//                adjacent_points[(k+2)%3] = Vs.row(index[(k+2)%3]);
//                real_pts[(k+2)%3] = true;
//            }
//            else
//            {
//                index[(k+2)%3] = -1; // Boundary edge
//                adjacent_points[(k+2)%3] << -1,-1,-1;
//                real_pts[(k+2)%3] = false;
//            }
//        }
//
//        GeoFeature::calculate_second_fundamental_form(Vs.row(Fs(i,0)), Vs.row(Fs(i,1)), Vs.row(Fs(i,2)), adjacent_points[0], adjacent_points[1], adjacent_points[2], real_pts, II);
//
//        II_list.push_back(I.inverse()*II);
//        std::cout<<"#"<<i<<std::endl;
//        std::cout<<I.inverse()*II<<std::endl;
//
//    }
//
//}
//
//void test_plane()
//{
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//
//    Eigen::MatrixXd V0;
//    Eigen::MatrixXi F0;
//    //Load a mesh in OFF format
//
//    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/rect-coarse.obj", V0, F0);
//    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/test_1.obj", V, F);
//    double _YoungsModulus = 1e4;
//    double _PossionRatio = 0.3;
//    double _thickness = 1;
//    auto op_shell = std::make_unique<ShellEnergy>();
//    double E_streching(0), E_bending(0);
//    Eigen::VectorXd diff_f_streching,diff_f_gravity,diff_f_bending;
//    op_shell->streching_energy(V,V0,F,_YoungsModulus,_PossionRatio,_thickness,E_streching,diff_f_streching);
//    op_shell->bending_energy(V, V0, F, _YoungsModulus, _PossionRatio, _thickness, E_bending, diff_f_bending);
//    double f = E_streching + E_bending;
//    std::cout<<f<<" "<<E_streching<<" "<<E_bending<<std::endl<<std::endl;
//    std::cout<<diff_f_streching.norm()<<std::endl<<std::endl;
//    std::cout<<diff_f_bending.norm()<<std::endl<<std::endl;
//    std::cout<<(diff_f_bending+diff_f_streching).norm()<<std::endl<<std::endl;
//}
//

#include <igl/opengl/glfw/Viewer.h>
#include "lib/ShellEnergy.h"
#include "lib/ShellSimulation.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>


auto op = std::make_unique<ShellSimulation>();
Eigen::MatrixXd V;
Eigen::MatrixXi F;


void reset()
{
    std::cout << std::endl << "Reset" << std::endl << std::endl;
    V = op->VU;
    F = op->F;
    
}

void repaint(igl::opengl::glfw::Viewer &viewer)
{
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        Eigen::MatrixXd colors(op->VU.rows(), 3);
        colors.col(0).setConstant(1.0);
        colors.col(1).setConstant(1.0);
        colors.col(2).setConstant(0.0);
        for (int i=0;i<op->p_fixed_index.size();i++)
        {
            colors(op->p_fixed_index[i], 0) = 1.0;
            colors(op->p_fixed_index[i], 1) = 0.0;
            colors(op->p_fixed_index[i], 2) = 0.0;
        }
        viewer.data().set_colors(colors);
}

int main(int argc, char *argv[])
{
    std::string problem = "/Users/chenzhen/UT/Research/shellbenchmarks/benchmarks/DrapedRect/3876_triangles/draped_rect";
    
        bool ok = op->set_up_simulation(problem);
        if (!ok)
        {
            std::cerr << "Couldn't load problem: " << problem << std::endl;
            return -1;
        }
    
    reset();
    
    igl::opengl::glfw::Viewer viewer;
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    // Add content to the default menu window
    menu.draw_viewer_menu_func = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();
        
        // Add new group
        if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset", ImVec2(-1, 0)))
            {
                reset();
                repaint(viewer);
            }
            if (ImGui::Button("Optimize Some Step", ImVec2(-1,0)))
            {
                op->compute_deformed_surface(V,F);
            }
        }
    };
    
    
    viewer.data().set_face_based(false);
    repaint(viewer);
    viewer.launch();
}

