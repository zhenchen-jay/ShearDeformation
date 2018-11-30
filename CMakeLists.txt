cmake_minimum_required(VERSION 3.3)
project(ShearDeformation)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

#########################################################################
# Use C++14
set (CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)

set(LIBIGL_INCLUDE_DIR ~/build/libigl-latest/include/)
#########################################################################
# libigl
# Add include files for libigl and Eigen
include_directories( ${LIBIGL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}
	/home/zxy/build/libigl-latest/external/glad/include/
	/home/zxy/build/libigl-latest/external/imgui/  )

# Self-defined libs
add_subdirectory(lib)
add_subdirectory(external/alglib)
add_subdirectory(external/LSQR)

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES} viewer3d.cpp)
target_link_libraries(${PROJECT_NAME}_bin lib alglib lsqr)

################################################################################
