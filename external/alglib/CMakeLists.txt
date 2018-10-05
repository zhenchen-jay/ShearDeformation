################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)

################################################################################

add_library(alglib
        ap.cpp
        ap.h
        alglibinternal.h
        alglibinternal.cpp
        alglibmisc.h
        alglibmisc.cpp
        dataanalysis.h
        dataanalysis.cpp
        diffequations.h
        diffequations.cpp
        fasttransforms.h
        fasttransforms.cpp
        integration.h
        integration.cpp
        interpolation.h
        interpolation.cpp
        linalg.h
        linalg.cpp
        optimization.h
        optimization.cpp
        solvers.h
        solvers.cpp
        specialfunctions.h
        specialfunctions.cpp
        statistics.h
        statistics.cpp
        stdafx.h)

set_target_properties(alglib PROPERTIES CXX_STANDARD 14)
set_target_properties(alglib PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(alglib PROPERTIES POSITION_INDEPENDENT_CODE ON)