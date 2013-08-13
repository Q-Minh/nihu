# cmake file for the NiHu project

# Set the minimum 
cmake_minimum_required (VERSION 2.6)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# Set the project name
project (NiHu)

set (NiHu_VERSION_MAJOR 1)
set (NiHu_VERSION_MINOR 0)

# Set the 3rdparty dir
set(NIHU_THIRDPARTY_DIR "${CMAKE_SOURCE_DIR}/ThirdParty")

# Find the number of system bits
include(Sysbits)

get_system_bits(NIHU_SYS_BITS)

message(STATUS "${NIHU_SYS_BITS}-bit system detected")

# COMPILER PROPERTIES
# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
	else ()
		message(STATUS "Compiler supports C++11 - OK")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()


### EIGEN section

# Eigen as headers is default in windows mode
if(WIN32)
	set(NIHU_EIGEN_AS_HEADERS 1)
endif(WIN32)

include(NiHuEigen)

### MATLAB Section

# look for environment variable MATLAB_ROOT
if (DEFINED ENV{MATLAB_ROOT})
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
elseif(DEFINED NIHU_MATLAB_ROOT)
	set(MATLAB_ROOT "${NIHU_MATLAB_ROOT}")
endif(DEFINED ENV{MATLAB_ROOT})

# look for matlabmex pacakge
find_package (MatlabMex REQUIRED)

# check if matlab mex is found
if(MATLABMEX_FOUND)
 	message(STATUS "Matlab MEX found, enabling build of mex files")
 	set(NIHU_BUILD_MEX 1)
	# force usage of matlab's mex comipler
	if(NIHU_FORCE_MEX_COMPILER)
		message(STATUS "Forcing Matlab's MEX compiler (${MATLAB_MEX}) for mex files")
	else(NIHU_FORCE_MEX_COMPILER)
		message(STATUS "Will use the ${CMAKE_CXX_COMPILER_ID} compiler (${CMAKE_CXX_COMPILER}) for mex files")
	endif(NIHU_FORCE_MEX_COMPILER)
# if matlab mex is not found
else(MATLABMEX_FOUND)
	message(STATUS "Matlab MEX not found, disabling build of mex files")
	set(NIHU_BUILD_MEX 0)
endif(MATLABMEX_FOUND)

# setup compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wunused -pedantic -O3 -std=c++11")

# setup global include directories
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${EIGEN_INCLUDE_DIRS}")

### Test section

enable_testing()
add_subdirectory(test)

# add tests if not disabled
# if(NOT NIHU_DISABLE_TESTS)
# 	add_subdirectory(test)
# else(NOT NIHU_DISABLE_TESTS)
# 	add_subdirectory(test EXCLUDE_FROM_ALL)
# endif(NOT NIHU_DISABLE_TESTS)

