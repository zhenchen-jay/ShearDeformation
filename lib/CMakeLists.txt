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

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

include_directories( ${LIBIGL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

################################################################################

add_library(lib
        ShellEnergy.h
        ShellEnergy.cpp
        ShellEnergyWithSwellRatio.h
        ShellEnergyWithSwellRatio.cpp
        ShellEnergyStandard.h
        ShellEnergyStandard.cpp
        ShellSimulation.h
        ShellSimulation.cpp
        ExternalEnergy.h
        ExternalEnergy.cpp
        GeometryFeature.h
        GeometryFeature.cpp
        MeshConnection.h
        MeshConnection.cpp
        elliptic_functions.cpp
        elliptic_functions.h
        ComputeCoefficient.cpp
        ComputeCoefficient.h
	    Equadprog.h
	    Equadprog.cpp
	    LSQR.h
	    LSQR.cpp
        FindFirstFundamentalCoef.h
        FindFirstFundamentalCoef.cpp)

set_target_properties(lib PROPERTIES CXX_STANDARD 14)
set_target_properties(lib PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(lib PUBLIC ..)

