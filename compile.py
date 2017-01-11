import os

COMPILE_HOKUSAI = "cmake ../"
COMPILE_EXAMPLE = "cmake ../"

os.system("cd hokusai/ && rm -r build/")
os.system("cd hokusai/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

os.system("cd examples/ && rm -r build/")
os.system("cd examples/ && mkdir build && cd build && "+COMPILE_HOKUSAI+" && "+"make")

