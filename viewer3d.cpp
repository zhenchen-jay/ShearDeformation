#include <igl/opengl/glfw/Viewer.h>
#include "lib/ShellEnergy.h"
#include "lib/ShellSimulation.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include "lib/ComputeCoefficient.h"
#include "lib/FindFirstFundamentalCoef.h"
#include <memory>

#ifndef IT_NUM
#define IT_NUM 10
#endif

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

void compute_hypar()
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/DrapedRect/3876_triangles/draped_rect_geometry.obj", Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        Vo(i,2) = 16*Vo(i,1) * Vo(i,0);
    }
    igl::writeOBJ("/Users/chenzhen/UT/Research/Results/hypar.obj", Vo, Fo);
}

void compute_cylinder()
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/DrapedRect/3876_triangles/draped_rect_geometry.obj", Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.025;
        double x = R*sin(Vo(i,0)/R);
        double z = R*cos(Vo(i,0)/R) - R;
        Vo(i,0) = x;
        Vo(i,2) = z;
    }
    igl::writeOBJ("/Users/chenzhen/UT/Research/Results/cylinder.obj", Vo, Fo);
}

int main(int argc, char *argv[])
{
    auto op_1 = std::make_unique<FindFirstFundamentalCoef>();
    op_1->test_func_grad();
    
//    //std::string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/CantileverPlate/1040_triangles/cantilever_plate";
//    //std:: string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/SlitAnnulus/2907_trianges/slit_annular_plate";
//    std::string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/DrapedRect/3876_triangles/draped_rect";
//    std::string tar_prob = "/Users/chenzhen/UT/Research/Results/cylinder";
//    //std::string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/ClampedCylindricalShell/2496_triangles/cantilever_cylindrical_shell";
//    bool ok = op->set_up_simulation(problem, tar_prob);
//    if (!ok)
//    {
//        std::cerr << "Couldn't load problem: " << problem << std::endl;
//        return -1;
//    }
//  //op->compute_deformed_surface(V, F);
//
//    reset();
//    igl::readOBJ(tar_prob + "/cylinder.obj", V, F);
//    // op->add_noise();
//    for(int i=0;i<IT_NUM;i++)
//    {
//        double ratio = (i+1)*1.0/IT_NUM;
//        op->ratio = ratio;
//        std::cout<<ratio<<std::endl;
//        op->compute_deformed_surface(V,F);
//    }
//    std::cout<<"Finished"<<std::endl;
//    igl::writeOBJ(tar_prob + "/Result.obj", V, F);
    
    
//    //op->compute_deformed_surface(V, F);
//    igl::opengl::glfw::Viewer viewer;
//
//    // Attach a menu plugin
//    igl::opengl::glfw::imgui::ImGuiMenu menu;
//    viewer.plugins.push_back(&menu);
//
//    // Add content to the default menu window
//    menu.draw_viewer_menu_func = [&]()
//    {
//        // Draw parent menu content
//        menu.draw_viewer_menu();
//
//        // Add new group
//        if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen))
//        {
//            if (ImGui::Button("Reset", ImVec2(-1, 0)))
//            {
//                reset();
//                repaint(viewer);
//            }
//            if (ImGui::Button("Optimize Some Step", ImVec2(-1,0)))
//            {
//                //igl::readOBJ("/Users/chenzhen/UT/Research/Results/cylinder.obj", V, F);
//                //op->add_noise();
//                for(int i=0;i<IT_NUM;i++)
//                {
//                    double ratio = (i+1)*1.0/IT_NUM;
//                    op->ratio = ratio;
//                    std::cout<<ratio<<std::endl;
//                    op->compute_deformed_surface(V,F);
//                }
//                repaint(viewer);
//                std::cout<<"Finished"<<std::endl;
//            }
//        }
//    };
//
//
//    viewer.data().set_face_based(false);
//    repaint(viewer);
//    viewer.launch();
}

