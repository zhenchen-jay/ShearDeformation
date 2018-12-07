cmake_minimum_required(VERSION 3.3)
project(ShearDeformation)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_PREFIX_PATH "$ENV{HOME}/.local")

#######################i#################################################
# Use C++14
set (CMAKE_CXX_STANDARD 14)

#########################################################################
# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)
###########################################################################


find_package(LIBIGL REQUIRED)
message("Found LibIGL at ${LIBIGL_INCLUDE_DIR} from $ENV(LIBIGL}")

find_package(ifopt REQUIRED)
message("Found IFOPT at ${ifopt_INCLUDE_DIRS} from $ENV(ifopt}")

# Add include files for libigl and Eigen
include_directories( ${LIBIGL_INCLUDE_DIR} )
include_directories( ${ifopt_INCLUDE_DIRS} )

# Self-defined libs
add_subdirectory(lib)
add_subdirectory(external/alglib)
add_subdirectory(external/LSQR)

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES} viewer3d.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui lib alglib lsqr ifopt::ifopt_ipopt)

################################################################################
