cmake_minimum_required(VERSION 4.2)
project(Solver1D_Euler_SOD_MKII)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Solver1D_Euler_SOD_MKII
        src/Cellule.cpp
        src/Cellule.h
        src/Domaine.cpp
        src/Domaine.h
        src/NumericalSchemes.cpp
        src/NumericalSchemes.h
        src/StateVector.h
        main.cpp)
