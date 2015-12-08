import os

DEIGEN_INCLUDE_DIRS="-DEIGEN_INCLUDE_DIRS=./../../../lib/eigen-3.2.6" 
DMAGNET_INCLUDE_DIRS="-DMAGNET_INCLUDE_DIRS=./../../../lib/magnet/include"
COMPILE_HOKUSAI = "cmake "+DEIGEN_INCLUDE_DIRS+" "+DMAGNET_INCLUDE_DIRS+" ../"
os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")
