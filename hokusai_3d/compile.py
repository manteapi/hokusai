import os

DALJABR_INCLUDE_DIRS="-DALJABR_INCLUDE_DIRS=~/Depot/github/aljabr" 

COMPILE_HOKUSAI = "cmake "+DALJABR_INCLUDE_DIRS+" ../"
COMPILE_EXAMPLE = "cmake "+DALJABR_INCLUDE_DIRS+" ../"

os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

os.system("cd examples/ && rm -r build/")
os.system("cd examples/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

