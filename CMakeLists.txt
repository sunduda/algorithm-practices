cmake_minimum_required(VERSION 3.13)  # CMake version check
project(tensor_toolkit)               # Create project "linear_algebra_library"
set(CMAKE_CXX_STANDARD 17)            # Enable c++17 standard

# Add main.cpp file of project root directory as source file
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(SOURCE_FILES ${PROJECT_SOURCE_DIR}/main.cpp)

# Add included directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add executable target with source files listed in SOURCE_FILES variable
add_library(tensor SHARED ${PROJECT_SOURCE_DIR}/src/tensor.cpp)
set(TENSOR_TOOLKIT_LIBRARY tensor ${PROJECT_SOURCE_DIR}/lib/libtensor.so)
# Add executable target with source files listed in SOURCE_FILES variable
#link_directories(src)
add_executable(tensor_toolkit_testapp ${SOURCE_FILES})

# Find and link libraries
target_link_libraries(tensor_toolkit_testapp ${TENSOR_TOOLKIT_LIBRARY} ${tensor_exceptions})
