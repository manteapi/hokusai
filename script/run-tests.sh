#!/bin/sh

cd hokusai_3D/test
mkdir build
cd build
cmake ..
make
make test
