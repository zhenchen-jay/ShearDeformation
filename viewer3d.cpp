//
// Created by Zhen Chen on 9/17/18.
//

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/boundary_loop.h>
#include <imgui/imgui.h>
//#include "lib/simulation.h"
#include "lib/ShellEnergy.h"
#include "lib/ShellSimulation.h"
#include "lib/GeometryFeature.h"
#include <iostream>
#include <string>
#include <fstream>

void write_boundary(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    std::vector<std::vector<int>> L;
    igl::boundary_loop(F, L);
    std::ofstream in;
    in.open("boundary.txt",std::ios::trunc);
    for(int i=0;i<L.size();i++)
    {
        in<<i<<std::endl; // The i-th boundary
        for(int j=0;j<L[i].size();j++)
        {
            in<<V(L[i][j],0)<<" "<<V(L[i][j],1)<<" "<<V(L[i][j],2)<<std::endl;
        }
    }
    in.close();
}

void test_grad()
{
    //auto op = std::make_unique<Simulation>();
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/plane.obj", V, F);
    //op->initialize(V, F, 0.3, 1e5, 1);
    for(int i=0;i<3*V.rows();i++)
    {
        if(i%3==2 || i == 7 || i==10 || i == 6 || i == 9)
            continue;
        std::cout<<std::endl;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(3*V.rows());
        v(i)=1;
        //for(int k = 3;k<10;k++)
        //    op->test_function(v, V, pow(10,-k));
    }
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::Matrix2d A;
    Eigen::Matrix2d B;
    
    Eigen::MatrixXd V0=V;
    Eigen::MatrixXi F0=F;
    
    // auto op = std::make_unique<Simulation>();
    auto op1 = std::make_unique<ShellSimulation>();

    //Load a mesh in OFF format
    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/rect-coarse.obj", V, F);
    //igl::readOBJ("/Users/chenzhen/UT/Research/Projects/ShearDeformation/Models/result.obj", V0, F0);
    
    //op->test(V, F, 0.3, 1e5, 1, 3);
    //op1->compute_deformed_surface(V,F,0.3,1e5,1,3);
    GeoFeature::test_face_normal();
    // write_boundary(V, F);
    //test_grad();
    
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    // Customize the menu
    float floatVariable = 0.1f; // Shared between two menus

    // Add content to the default menu window
    menu.draw_custom_window_func = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputFloat("float", &floatVariable, 0, 0, 3);

            // ... or using a custom callback
            static bool boolVariable = true;
            if (ImGui::Checkbox("bool", &boolVariable))
            {
                // do something
                std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
            }

            // Expose an enumeration type
            enum Orientation { Up=0, Down, Left, Right };
            static Orientation dir = Up;
            ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

            // We can also use a std::vector<std::string> defined dynamically
            static int num_choices = 3;
            static std::vector<std::string> choices;
            static int idx_choice = 0;
            if (ImGui::InputInt("Num letters", &num_choices))
            {
                num_choices = std::max(1, std::min(26, num_choices));
            }
            if (num_choices != (int) choices.size())
            {
                choices.resize(num_choices);
                for (int i = 0; i < num_choices; ++i)
                    choices[i] = std::string(1, 'A' + i);
                if (idx_choice >= num_choices)
                    idx_choice = num_choices - 1;
            }
            ImGui::Combo("Letter", &idx_choice, choices);

            // Add a button
            if (ImGui::Button("Print Hello", ImVec2(-1,0)))
            {
                std::cout << "Hello\n";
            }
        }
    };

    // Draw additional windows
    menu.draw_custom_window_func = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
                "New Window", nullptr,
                ImGuiWindowFlags_NoSavedSettings
        );

        // Expose the same variable directly ...
        ImGui::PushItemWidth(-80);
        ImGui::DragFloat("float", &floatVariable, 0.0, 0.0, 3.0);
        ImGui::PopItemWidth();

        static std::string str = "bunny";
        ImGui::InputText("Name", str);

        ImGui::End();
    };
    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.launch();
}



