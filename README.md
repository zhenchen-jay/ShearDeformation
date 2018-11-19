# Shear Deformation project

This project is used for compute the deformed surface.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `ShearDeformation__bin` binary.

## Run

From within the `build` directory just issue:

    ./ShearDeformation_bin


## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

We recommend you to install libigl using git via:

    git clone https://github.com/libigl/libigl.git
    cd libigl/
    git checkout 6ebc585611d27d8e2038bbbf9cb4910c51713848
    git submodule update --init --recursive
    cd ..

If you have installed libigl at `/path/to/libigl/` then a good place to clone
this library is `/path/to/yourproject/`.
