################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
# project(lib)

################################################################################
# Use C++14
# Generate position independent code

################################################################################
# Required libraries
################################################################################

find_package(LIBIGL REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories( ${LIBIGL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})

################################################################################

add_library(lib
        simulation.h
        simulation.cpp
        ShellEnergy.h
        ShellEnergy.cpp
        ShellSimulation.h
        ShellSimulation.cpp
        GravityEnergy.h
        GravityEnergy.cpp
        GeometryFeature.h
        GeometryFeature.cpp
        MeshConnection.h
        MeshConnection.cpp
        )

set_target_properties(lib PROPERTIES CXX_STANDARD 14)
set_target_properties(lib PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(lib PUBLIC ..)

