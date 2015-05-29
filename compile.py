import os

DALJABR_INCLUDE_DIRS="-DALJABR_INCLUDE_DIRS=~/Depot/github/aljabr" 
DEIGEN_INCLUDE_DIRS="-DEIGEN_INCLUDE_DIRS=~/Depot/github/lib/eigen-3.2.4" 
DMAGNET_INCLUDE_DIRS="-DMAGNET_INCLUDE_DIRS=~/Depot/github/hokusai/hokusai/lib/magnet/include"

COMPILE_HOKUSAI = "cmake "+DALJABR_INCLUDE_DIRS+" "+DEIGEN_INCLUDE_DIRS+" "+DMAGNET_INCLUDE_DIRS+" ../"

os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

