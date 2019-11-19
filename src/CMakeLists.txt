# CMakeLists.txt
# cmake file for the NiHu project
# 
# This file and its submodules follow the coding syntax defined here:
#   https://community.kde.org/Policies/CMake_Coding_Style

# NOTE: There is a problem of Win10 + Matlab R2018a / symbolic math toolbox

###
### cmake and project section
###

# Set the minimum required cmake version
cmake_minimum_required (VERSION 3.15.0 FATAL_ERROR)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# Set the project name
project(NiHu)

set(NIHU_VERSION_MAJOR 2)
set(NIHU_VERSION_MINOR 0)

# Find the number of system bits
include(NiHuSystem)
# Configure compiler and linker
include(NiHuCompilerSettings)
include(NiHuLinkerSettings)

###
### ThirdParty software 
###

# Set the 3rdparty dir
set(NIHU_THIRDPARTY_DIR "${CMAKE_BINARY_DIR}/ThirdParty")

### Boost section
include(NiHuBoost)

### EIGEN section
include(NiHuEigen)

### FFTW section 
include(NiHuFFTW)

### MATLAB section
if (NOT DEFINED NIHU_DISABLE_MATLAB)
	include(NiHuMatlab)
endif()

### GTest section
if(NIHU_ENABLE_TESTING)
	# Check for gtest module
	include (NiHuGTest)
endif()
	
###
### Modules configuration 
###

### FMM module configuration
if(NOT DEFINED NIHU_DISABLE_FMM)
	include(NiHuFMM)
endif()

###
### Include setup for compilation
###

# setup global include directories
include_directories("${CMAKE_SOURCE_DIR}")
if(DEFINED Boost_INCLUDE_DIRS)
    include_directories("${Boost_INCLUDE_DIRS}")
endif()
include_directories("${EIGEN_INCLUDE_DIRS}")
include_directories("${FFTW3_INCLUDE_DIRS}")

###
### Global library targets
###
include(NiHuLibraries)

### Subdirectories
# Tutorials directory
if(NOT DEFINED NIHU_DISABLE_TUTORIALS)
	add_subdirectory(tutorial)
endif()

if(NOT DEFINED NIHU_DISABLE_APPLICATIONS)
	add_subdirectory(applications/helmholtz)
	add_subdirectory(applications/laplace)
endif()

### Test section
if(NIHU_ENABLE_TESTING)
	include(NiHuTesting)
endif()

### Documentation section
if(NOT DEFINED NIHU_DISABLE_DOC)
	include(NiHuDoc)
endif()

### Installation section
if(NOT DEFINED NIHU_DISABLE_INSTALL)
	include(NiHuInstall)
endif()
