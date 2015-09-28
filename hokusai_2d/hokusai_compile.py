import os

DEIGEN_INCLUDE_DIRS="-DEIGEN_INCLUDE_DIRS=./../../extlib/eigen-3.2.5" 
DMAGNET_INCLUDE_DIRS="-DMAGNET_INCLUDE_DIRS=./../../extlib/magnet/include"
COMPILE_HOKUSAI = "cmake "+DEIGEN_INCLUDE_DIRS+" "+DMAGNET_INCLUDE_DIRS+" ../"
os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")
