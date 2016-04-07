#!/bin/sh

cd hokusai_3d/test
mkdir build
cd build
cmake ..
make
make test
