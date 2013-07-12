# cmake file for the NiHu project

# Set the minimum 
cmake_minimum_required (VERSION 2.6)

project (NiHu)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# Check for eigen
find_package (Eigen REQUIRED)

# Check if eigen found
if(EIGEN_FOUND)
	include_directories ("${EIGEN_INCLUDE_DIRS}")
else(EIGEN_FOUND)
	message(FATAL_ERROR "Eigen3 not found, CMake will exit." )
endif(EIGEN_FOUND)

# look for environment variable MATLAB_ROOT
if (DEFINED ENV{MATLAB_ROOT})
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
endif(DEFINED ENV{MATLAB_ROOT})

# look for matlabmex pacakge
find_package (MatlabMex REQUIRED)

if(MATLABMEX_FOUND)
 	message(STATUS "Matlab MEX found, enabling buid of mex files")
 	set(BUILD_MEX 1)
else(MATLABMEX_FOUND)
	MESSAGE(STATUS "Matlab MEX not found, disabling build of mex files")
	set(BUILD_MEX 0)
endif(MATLABMEX_FOUND)

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -pedantic -O3 -std=c++11" )

add_subdirectory(test)
