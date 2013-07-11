# cmake file for the NiHu project

# Set the minimum 
cmake_minimum_required (VERSION 2.6)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# Check for eigen
find_package (Eigen REQUIRED)

# Check if eigen found
if(EIGEN_FOUND)
	include_directories ("${EIGEN_INCLUDE_DIRS}")
else()
	message(FATAL_ERROR "Eigen3 not found, CMake will exit." )
endif()

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -std=c++11" )

project (NiHu)
add_subdirectory(test)
# add_executable(Tutorial tutorial.cxx)