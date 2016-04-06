import os

DALJABR_INCLUDE_DIRS="-DALJABR_INCLUDE_DIRS=~/Depot/github/aljabr" 
DEIGEN_INCLUDE_DIRS="-DEIGEN_INCLUDE_DIRS=./../../extlib/eigen-3.2.5" 
DMAGNET_INCLUDE_DIRS="-DMAGNET_INCLUDE_DIRS=./../../extlib/magnet/include"
COMPILE_HOKUSAI = "cmake "+ DALJABR_INCLUDE_DIRS+" "+DEIGEN_INCLUDE_DIRS+" "+DMAGNET_INCLUDE_DIRS+" ../"
#os.system("cd extlib/AntTweakBar/src/ && make")
#os.system("cd extlib/glew-1.10.0/src/ && make")
#os.system("cd extlib/sfml-2.1 && mkdir build && cd build && cmake .. && make")
os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")
