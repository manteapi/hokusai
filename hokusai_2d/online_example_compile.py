import os

DEIGEN_INCLUDE_DIRS="-DEIGEN_INCLUDE_DIRS=./../../extlib/eigen-3.2.5" 
DMAGNET_INCLUDE_DIRS="-DMAGNET_INCLUDE_DIRS=./../../extlib/magnet/include"
GLM_INCLUDE_DIRS="-DGLM_INCLUDE_DIRS=./../../extlib/glm-0.9.4.0"
ANTTWEAKBAR_INCLUDE_DIRS="-DANTTWEAKBAR_INCLUDE_DIRS=./../../extlib/AntTweakBar/include"
ANTTWEAKBAR_LIBRARIES="-DANTTWEAKBAR_LIBRARIES=./../../extlib/AntTweakBar/lib/libAntTweakBar.so"
GLEW_INCLUDE_DIRS="-DGLEW_INCLUDE_DIRS=./../../extlib/glew-1.10.0"
GLEW_LIBRARIES="-DGLEW_LIBRARIES=./../../extlib/glew-1.10.0/lib/libglew.so"
SFML_INCLUDE_DIRS="-DSFML_INCLUDE_DIRS=./../../extlib/sfml-2.1/include"
SFML_LIBRARIES_GRAPHICS="-DSFML_GRAPHICS_LIBRARIES=./../../extlib/sfml-2.1/build/lib/libsfml-graphics.so"
SFML_LIBRARIES_WINDOW="-DSFML_WINDOW_LIBRARIES=./../../extlib/sfml-2.1/build/lib/libsfml-window.so"
SFML_LIBRARIES_SYSTEM="-DSFML_SYSTEM_LIBRARIES=./../../extlib/sfml-2.1/build/lib/libsfml-system.so"
HOKUSAI_OPENMP="-DHOKUSAI_OPENMP=true"
HOKUSAI_INCLUDE_DIRS="-DHOKUSAI_INCLUDE_DIRS=./../../hokusai/include"
HOKUSAI_LIBRARIES="-DHOKUSAI_LIBRARIES=./../../hokusai/build/libHOKUSAI.a"

COMPILE_ONLINE_EXAMPLE = "cmake "+DEIGEN_INCLUDE_DIRS+" "+DMAGNET_INCLUDE_DIRS+" "+GLM_INCLUDE_DIRS+" "+GLEW_INCLUDE_DIRS+" "+GLEW_LIBRARIES+" "+SFML_INCLUDE_DIRS+" "+SFML_LIBRARIES_GRAPHICS+" "+SFML_LIBRARIES_WINDOW+" "+SFML_LIBRARIES_SYSTEM+" "+ANTTWEAKBAR_INCLUDE_DIRS+" "+ANTTWEAKBAR_LIBRARIES+" "+HOKUSAI_INCLUDE_DIRS+" "+HOKUSAI_LIBRARIES+" "+HOKUSAI_OPENMP+" ../"

os.system("cd online_example/ && rm -r build/")
os.system("cd online_example/ && mkdir build && cd build && "+COMPILE_ONLINE_EXAMPLE+" && "+"make")

