HOKUSAI - Liquid simulation using SPH
=====================================

A C++ implementation of the implicit incompressible sph model (IISPH) proposed by Ihmsen et al. 2014.

The original paper can be found here : http://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf

How to compile the library
----------------------------
    cd hokusai
    mkdir build
    cmake ..
    make

How to compile/run examples
------------------------------------
    cd examples
    mkdir build
    cd build
    cmake ..
    make
    ./dambreak

How to compile/run tests
--------------------------
    cd tests
    mkdir build
    cd build
    cmake ..
    make test

How to generate the documentation
-----------------------------------
    cd hokusai
    mkdir build
    cmake ..
    make doc
