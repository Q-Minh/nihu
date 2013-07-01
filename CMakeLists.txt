# cmake file for the NiHu project

# Set the minimum 
cmake_minimum_required (VERSION 2.6)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Check for eigen
find_package (Eigen REQUIRED)


project (NiHu)
# add_executable(Tutorial tutorial.cxx)