cmake_minimum_required(VERSION 3.10)

# Project Name
project(LFEM_Quad4)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/eigen ${CMAKE_SOURCE_DIR}/include)

# Add source files
set(SOURCE_FILES
    src/main.cpp
    src/Node.cpp
    src/Element.cpp
    src/FEMSolver.cpp
)

# Create executable
add_executable(LFEM_Quad4 ${SOURCE_FILES})