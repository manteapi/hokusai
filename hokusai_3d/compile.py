import os

DALJABR_INCLUDE_DIRS="-DALJABR_INCLUDE_DIRS=~/Depot/github/aljabr" 

DHOKUSAI_INCLUDE_DIRS="-DHOKUSAI_INCLUDE_DIRS=~/Depot/github/hokusai/hokusai_3d/hokusai/include/" 
DHOKUSAI_LIBRARIES="-DHOKUSAI_LIBRARIES=~/Depot/github/hokusai/hokusai_3d/hokusai/build/libHOKUSAI.a" 

COMPILE_HOKUSAI = "cmake "+DALJABR_INCLUDE_DIRS+" ../"
COMPILE_EXAMPLE = "cmake "+DALJABR_INCLUDE_DIRS+" "+DHOKUSAI_INCLUDE_DIRS+" "+DHOKUSAI_LIBRARIES+" ../"

os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

os.system("cd examples/ && rm -r build/")
os.system("cd examples/ && mkdir build && cd build && "+COMPILE_EXAMPLE+" && "+"make")

