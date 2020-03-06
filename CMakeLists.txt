cmake_minimum_required(VERSION 3.13)  # CMake version check
project(tensor_toolkit)               # Create project "linear_algebra_library"
set(CMAKE_CXX_STANDARD 17)            # Enable c++17 standard

# Add main.cpp file of project root directory as source file
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib )
set( SOURCE_FILES main.cpp)

# Add included directories
include_directories(includes)

# Add executable target with source files listed in SOURCE_FILES variable
add_library( tensor_toolkit SHARED src/tensor.cpp)

# Add executable target with source files listed in SOURCE_FILES variable
#link_directories(src)
add_executable( tensor_toolkit_testapp ${SOURCE_FILES})

# Find and link libraries
#find_library(TENSOR_TOOLKIT tensor_toolkit ${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(tensor_toolkit_testapp ${tensor_toolkit})
