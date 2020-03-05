cmake_minimum_required(VERSION 3.13)  # CMake version check
project(tensor_toolkit)               # Create project "linear_algebra_library"
set(CMAKE_CXX_STANDARD 17)            # Enable c++17 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp include/xinda_math/tensor.cpp)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(tensor_toolkit ${SOURCE_FILES})