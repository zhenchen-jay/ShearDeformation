#include <memory>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/barycenter.h>
#include <imgui/imgui.h>
#include <igl/opengl/MeshGL.h>
#include <Eigen/Eigenvalues>
// #include <ifopt/problem.h>
// #include <ifopt/ipopt_solver.h>
// #include <ifopt/test_vars_constr_cost.h>

#include "lib/ComputeCoefficient.h"
#include "lib/FindFirstFundamentalCoef.h"
#include "lib/ShellEnergy.h"
#include "lib/ShellSimulation.h"
#include "lib/GeometryFeature.h"


// #include "lib/IfoptSolver.h"

#ifndef IT_NUM
#define IT_NUM 1
#endif

auto op = std::make_unique<ShellSimulation>();
Eigen::MatrixXd V;
Eigen::MatrixXi F;
std::vector<double> L_list;

bool _is_show_vec = false;
bool _is_show_f_colors = false;

void load_L_list(std::string file_path)
{
    std::ifstream infile(file_path);
    if(!infile)
        return;
    int num;
    infile >> num;
    double d;
    for(int i=0;i<num;i++)
    {
        infile >> d;
        L_list.push_back(d);
    }
}

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
    colors.col(2).setConstant(0);
//    for (int i=0;i<op->p_fixed_index.size();i++)
//    {
//        colors(op->p_fixed_index[i], 0) = 1.0;
//        colors(op->p_fixed_index[i], 1) = 0.0;
//        colors(op->p_fixed_index[i], 2) = 0.0;
//    }
    viewer.data().set_colors(colors);
    viewer.data().line_width = 2;
    viewer.core.background_color<<1,1,1,1;
    
    if(_is_show_vec)
    {
        Eigen::MatrixXd BC, Vec1, Vec2;
        igl::barycenter(V, F, BC);
        Vec1.resize(F.rows(), 3);
        Vec2.resize(F.rows(), 3);
        load_L_list("../../benchmarks/TestModels/L_list_cylinder.dat");
        for(int i=0;i<F.rows();i++)
        {
            Eigen::Matrix2d ID, IU, A, IM;
            GeoFeature::calculate_first_fundamental_form(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),ID);
            GeoFeature::calculate_first_fundamental_form(op->VU.row(F(i,0)),op->VU.row(F(i,1)),op->VU.row(F(i,2)),IU);
            IM<< L_list[3*i],0,
            L_list[3*i+1],L_list[3*i+2];
            IM = IM*IM.transpose();
            A = IU.inverse()*ID;
            //A = IU.inverse()*IM;
            std::cout<<A<<std::endl;
            Eigen::EigenSolver<Eigen::MatrixXd> es(A);
            double eigen_value_1 = es.eigenvalues()[0].real();
            Eigen::VectorXd eigen_vec_1 = es.eigenvectors().col(0).real();
            
            double eigen_value_2 = es.eigenvalues()[1].real();
            Eigen::VectorXd eigen_vec_2 = es.eigenvectors().col(1).real();
            
            std::cout<<eigen_value_1<<" "<<eigen_value_2<<std::endl;
            break;
            
            Vec1.row(i) = eigen_value_1*(eigen_vec_1(0)*(V.row(F(i,1))-V.row(F(i,0))) + eigen_vec_1(1)*(V.row(F(i,2))-V.row(F(i,0))));
            Vec2.row(i) = eigen_value_2*(eigen_vec_2(0)*(V.row(F(i,1))-V.row(F(i,0))) + eigen_vec_2(1)*(V.row(F(i,2))-V.row(F(i,0))));
            
        }
        const Eigen::RowVector3d red(1,0,0), black(0,0,0);
        viewer.data().add_edges(BC,BC+Vec1, red);
        viewer.data().add_edges(BC,BC+Vec2, black);
    }
    
    if(_is_show_f_colors)
    {
        load_L_list("../../benchmarks/TestModels/L_list_cylinder.dat");
        int faces = F.rows();
        igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_PARULA;
        Eigen::VectorXd Z(faces);
        Eigen::MatrixXd faceColors(faces, 3);
        
        for(int i=0;i<F.rows();i++)
        {
            Eigen::Matrix2d ID, IU, A, IM;
            GeoFeature::calculate_first_fundamental_form(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),ID);
            GeoFeature::calculate_first_fundamental_form(op->VU.row(F(i,0)),op->VU.row(F(i,1)),op->VU.row(F(i,2)),IU);
            IM<< L_list[3*i],0,
            L_list[3*i+1],L_list[3*i+2];
            IM = IM*IM.transpose();
            A = IU.inverse()*ID;
            //A = IU.inverse()*IM;
            Eigen::EigenSolver<Eigen::MatrixXd> es(A);
            double eigen_value_1 = es.eigenvalues()[0].real();
            
            double eigen_value_2 = es.eigenvalues()[1].real();
            Z(i) = eigen_value_1 + eigen_value_2;
            std::cout<<eigen_value_1<<" "<<eigen_value_2<<std::endl;
            break;
            
        }
        
        igl::colormap(viz_color, Z, true, faceColors); // true here means libigl will automatically normalize Z, which may or may not be what you want.
        viewer.data().set_colors(faceColors);
    }
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

void compute_sphere()
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ("../../benchmarks/DrapedRect/3876_triangles/draped_rect_geometry.obj", Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.01;
        double u = Vo(i,0);
        double v = Vo(i,1);
        double z = R - R*R/sqrt(R*R+u*u+v*v);
        double x = (R-z)/R*u;
        double y = (R-z)/R*v;
        Vo(i,0) = x;
        Vo(i,1) = y;
        Vo(i,2) = z;
    }
    igl::writeOBJ("../../benchmarks/sphere.obj", Vo, Fo);
    
}


int main(int argc, char *argv[])
{
//    compute_sphere();
//    return 0;
//    using namespace ifopt;
//    double YoungsModulus = 1e5;
//    double PossionRatio = 0.3;
//    double thickness = 1e-3;
//
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//    igl::readOBJ("../../benchmarks/TestModels/sphere.obj", V, F);
//    Problem nlp;
//    auto variable_set = std::make_shared<optVariables>(3*F.rows(), "var_set");
//    auto constraint_set = std::make_shared<optConstraint>(F.rows(), "constraint");
//    auto cost_term = std::make_shared<optCost>(V,F,YoungsModulus,PossionRatio,thickness);
//    nlp.AddVariableSet(variable_set);
//    nlp.AddConstraintSet(constraint_set);
//    nlp.AddCostSet(cost_term);
//    nlp.PrintCurrent();
//
//    IpoptSolver ipopt;
//    ipopt.SetOption("linear_solver", "mumps");
//    ipopt.SetOption("jacobian_approximation", "exact");
//    ipopt.SetOption("max_cpu_time", 1e6);
//    ipopt.SetOption("tol", 1e-10);
//    ipopt.SetOption("print_level", 5);
//
//    // 3 . solve
//    ipopt.Solve(nlp);
//    Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
//    //std::cout << x.transpose() << std::endl;
//    auto op = std::make_unique<FindFirstFundamentalCoef>();
//    op->set_up(V, F, YoungsModulus, PossionRatio, thickness);
//    double f;
//    Eigen::VectorXd df;
//    op->get_func_grad(x, f, df);
//    std::cout<<f<<std::endl;
//    std::cout<<df<<std::endl;

    
    
//    auto op_1 = std::make_unique<FindFirstFundamentalCoef>();
//    op_1->test_func_grad();
    
    //std::string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/CantileverPlate/1040_triangles/cantilever_plate";
    //std:: string problem = "/Users/chenzhen/UT/Research/Projects/ShearDeformation/benchmarks/SlitAnnulus/2907_trianges/slit_annular_plate";
    std::string problem = "../../benchmarks/DrapedRect/3876_triangles/draped_rect";
    std::string tar_prob = "../../benchmarks/TestModels";
    bool ok = op->set_up_simulation(problem, tar_prob);
    if (!ok)
    {
        std::cerr << "Couldn't load problem: " << problem << std::endl;
        return -1;
    }
  //op->compute_deformed_surface(V, F);

    reset();
    igl::readOBJ(tar_prob + "/cylinder.obj",V,F);
    //tar_prob + "/sphere.obj""", V, F);
    // op->add_noise();
//    for(int i=0;i<IT_NUM;i++)
//    {
//        double ratio = (i+1)*1.0/IT_NUM;
//        op->ratio = ratio;
//        std::cout<<ratio<<std::endl;
//        op->compute_deformed_surface(V,F);
//    }
//    std::cout<<"Finished"<<std::endl;
//    igl::writeOBJ(tar_prob + "/Result.obj", V, F);
    
    
    //op->compute_deformed_surface(V, F);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    igl::opengl::glfw::Viewer viewer;
    viewer.plugins.push_back(&menu);

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Workspace
        if (ImGui::CollapsingHeader("Workspace", ImGuiTreeNodeFlags_DefaultOpen))
        {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            if (ImGui::Button("Load##Workspace", ImVec2((w-p)/2.f, 0)))
            {
                viewer.load_scene();
            }
            ImGui::SameLine(0, p);
            if (ImGui::Button("Save##Workspace", ImVec2((w-p)/2.f, 0)))
            {
                viewer.save_scene();
            }
        }
        
        // Mesh
        if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
        {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            if (ImGui::Button("Load##Mesh", ImVec2((w-p)/2.f, 0)))
            {
                viewer.data().clear();
                viewer.open_dialog_load_mesh();
                V = viewer.data().V;
                F = viewer.data().F;
                repaint(viewer);
            }
            ImGui::SameLine(0, p);
            if (ImGui::Button("Save##Mesh", ImVec2((w-p)/2.f, 0)))
            {
                viewer.open_dialog_save_mesh();
            }
        }
        
        if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Center object", ImVec2(-1, 0)))
            {
                viewer.core.align_camera_center(viewer.data().V, viewer.data().F);
            }
            if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
            {
                viewer.snap_to_canonical_quaternion();
            }
            
            // Select rotation type
            int rotation_type = static_cast<int>(viewer.core.rotation_type);
            static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
            static bool orthographic = true;
            if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
            {
                using RT = igl::opengl::ViewerCore::RotationType;
                auto new_type = static_cast<RT>(rotation_type);
                if (new_type != viewer.core.rotation_type)
                {
                    if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
                    {
                        trackball_angle = viewer.core.trackball_angle;
                        orthographic = viewer.core.orthographic;
                        viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
                        viewer.core.orthographic = true;
                    }
                    else if (viewer.core.rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
                    {
                        viewer.core.trackball_angle = trackball_angle;
                        viewer.core.orthographic = orthographic;
                    }
                    viewer.core.set_rotation_type(new_type);
                }
            }
            
            // Orthographic view
            ImGui::Checkbox("Orthographic view", &(viewer.core.orthographic));
            ImGui::PopItemWidth();
        }
        
        // Draw options
        if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if(ImGui::Checkbox("Vector Field", &_is_show_vec))
            {
                repaint(viewer);
            }
            if(ImGui::Checkbox("Face Colors", &_is_show_f_colors))
            {
                repaint(viewer);
            }
            ImGui::Checkbox("Show texture", &(viewer.data().show_texture));
            if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals)))
            {
                viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
            }
            ImGui::Checkbox("Show overlay", &(viewer.data().show_overlay));
            ImGui::Checkbox("Show overlay depth", &(viewer.data().show_overlay_depth));
            ImGui::ColorEdit4("Background", viewer.core.background_color.data(),
                              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
            ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
                              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
            ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
            ImGui::PopItemWidth();
        }
        
        // Overlays
        if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
            ImGui::Checkbox("Fill", &(viewer.data().show_faces));
            ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
            ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
        }

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
//                igl::readOBJ("../../benchmarks/TestModels/cylinder.obj", V, F);
                //op->add_noise();
                for(int i=0;i<IT_NUM;i++)
                {
                    double ratio = (i+1)*1.0/IT_NUM;
                    op->ratio = ratio;
                    std::cout<<ratio<<std::endl;
                    op->compute_deformed_surface(V,F);
                }
                repaint(viewer);
                std::cout<<"Finished"<<std::endl;
            }
        }
    };


    viewer.data().set_face_based(false);
    repaint(viewer);
    viewer.launch();
 }

