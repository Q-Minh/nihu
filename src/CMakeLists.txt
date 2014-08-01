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
	message(STATUS "Setting up additional flags for Clang")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

	
# endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

### Linker configuration
if(WIN32)
	SET(CMAKE_EXE_LINKER_FLAGS "-MD ${CMAKE_EXE_LINKER_FLAGS} -s")
	SET(CMAKE_SHARED_LINKER_FLAGS "-MD ${CMAKE_SHARED_LINKER_FLAGS} -s")
endif(WIN32)

### EIGEN section

# Perform eigen installation
include(NiHuEigen)

### MATLAB Section

# look for environment variable MATLAB_ROOT
if (DEFINED ENV{MATLAB_ROOT})
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
elseif(DEFINED NIHU_MATLAB_PATH)
	set(MATLAB_ROOT "${NIHU_MATLAB_PATH}")
endif(DEFINED ENV{MATLAB_ROOT})

# look for matlabmex pacakge
find_package (MatlabMex)

# check if matlab mex is found
if(MATLABMEX_FOUND)
 	message(STATUS "Matlab MEX found, enabling build of mex files")
 	set(NIHU_BUILD_MEX 1)
	# force usage of matlab's mex comipler
	if(NIHU_MATLAB_FORCE_MEX_COMPILER)
		message(STATUS "Forcing Matlab's MEX compiler (${MATLAB_MEX}) for mex files")
	else(NIHU_MATLAB_FORCE_MEX_COMPILER)
		message(STATUS "Will use the ${CMAKE_CXX_COMPILER_ID} compiler (${CMAKE_CXX_COMPILER}) for mex files")
	endif(NIHU_MATLAB_FORCE_MEX_COMPILER)
	# set the flags
	if(WIN32)
		SET(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		SET(MEX_SHARED_LINKER_FLAGS "-shared ${CMAKE_SHARED_LINKER_FLAGS} -m${NIHU_SYS_BITS} -L\"${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}\" -L\"${MATLAB_ROOT}/extern/lib/win${NIHU_SYS_BITS}/microsoft\" -L\"${CMAKE_BINARY_DIR}/lib\"  -lmex -lmx -lmat")
	else(WIN32)
		set(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DMATLAB_MEX_FILE")
	endif(WIN32)
	# Set the extension for mex files
	set(MEX_SHARED_LIBRARY_SUFFIX "${MATLAB_MEXEXT}")
	# Add the include directory
	include_directories("${MATLAB_ROOT}/extern/include")
# if matlab mex is not found
else(MATLABMEX_FOUND)
	message(STATUS "Matlab MEX not found, disabling build of mex files")
	set(NIHU_BUILD_MEX 0)
endif(MATLABMEX_FOUND)

### Setup directories for installation
if(NOT DEFINED NIHU_INSTALL_DIR)
	set(NIHU_INSTALL_DIR ${CMAKE_BINARY_DIR})
endif(NOT DEFINED NIHU_INSTALL_DIR)

set(CMAKE_INSTALL_PREFIX ${NIHU_INSTALL_DIR})

# setup compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wunused -pedantic -O3 -std=c++11")

# setup global include directories
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${EIGEN_INCLUDE_DIRS}")

### Libraries
if(WIN32)
	set(NIHU_STATIC_LIB_SUFFIX ".lib")
else()
	set(NIHU_STATIC_LIB_SUFFIX ".a")
endif()
	
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

FILE(GLOB NIHU_COMMON_LIBRARIES 
	"${CMAKE_SOURCE_DIR}/library/lib_*.cpp")

# Iterate through the libraries by foreach
foreach(lib_source ${NIHU_COMMON_LIBRARIES})
	# Extract relative path
	file(RELATIVE_PATH local_lib_source "${CMAKE_SOURCE_DIR}/library" ${lib_source})
	string(REPLACE ".cpp" "" lib_name ${local_lib_source})
	# Add the library
	add_library(${lib_name} STATIC ${lib_source})
	# Set the target properties
	set_target_properties(${lib_name}
		PROPERTIES 
		PREFIX "" 
		SUFFIX ${NIHU_STATIC_LIB_SUFFIX}
		LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	)
	# Set the link libraries
	SET(NIHU_LINK_LIBRARIES ${NIHU_LINK_LIBRARIES}; ${lib_name})
endforeach(lib_source)

### Subdirectories
add_subdirectory(tutorial)

### Test section
if(NIHU_ENABLE_TESTING)
	if(NIHU_ENABLE_TEST_INSTALL)
		if(${NIHU_INSTALL_DIR} MATCHES ${CMAKE_BINARY_DIR})
			set(NIHU_ENABLE_TEST_INSTALL 0)
			message(STATUS "Installation of tests disabled (build dir is same as install dir)")
		endif()
	endif(NIHU_ENABLE_TEST_INSTALL)
	enable_testing()
	if(NOT NIHU_DISABLE_TEST_BUILD)
		add_subdirectory(test)
	else(NOT NIHU_DISABLE_TEST_BUILD)
		add_subdirectory(test EXCLUDE_FROM_ALL)
	endif(NOT NIHU_DISABLE_TEST_BUILD)
endif(NIHU_ENABLE_TESTING)

### Documentation section

find_package(Doxygen)

if(DOXYGEN_FOUND)
	# Check documentation path
	if(DEFINED NIHU_HTML_DOC_DIR)
	else(DEFINED NIHU_HTML_DOC_DIR)
		set(NIHU_HTML_DOC_DIR "${CMAKE_BINARY_DIR}/doc/html")
		file(MAKE_DIRECTORY ${NIHU_HTML_DOC_DIR})
	endif(DEFINED NIHU_HTML_DOC_DIR)
	message(STATUS "html documentation path: ${NIHU_HTML_DOC_DIR}")
	
	if(NOT NIHU_MATHJAX_DISABLE)
		message(STATUS "Enabling MathJax in Doxygen documentation")
		set(NIHU_MATHJAX_USE "YES")
		if(DEFINED NIHU_MATHJAX_PATH)
			file(RELATIVE_PATH NIHU_MATHJAX_RELPATH "${NIHU_HTML_DOC_DIR}" "${NIHU_MATHJAX_PATH}")
			message(STATUS "MathJax relative path: ${NIHU_MATHJAX_RELPATH}")
		else(DEFINED NIHU_MATHJAX_PATH)
			set(NIHU_MATHJAX_RELPATH "http://cdn.mathjax.org/mathjax/latest")
			message(STATUS "Will use CDN for MathJax: ${NIHU_MATHJAX_RELPATH}")
		endif(DEFINED NIHU_MATHJAX_PATH)
	elseif(NOT NIHU_MATHJAX_DISABLE)
		message(STATUS "Disabling MathJax in Doxygen documentation")
		set(NIHU_MATHJAX_USE "NO")
		set(NIHU_MATHJAX_RELPATH "")
	endif(NOT NIHU_MATHJAX_DISABLE)

	# Configuration of doxyfiles
	configure_file(${CMAKE_SOURCE_DIR}/../Doxyfile.in 
		${CMAKE_BINARY_DIR}/Doxyfile)

	# doc target - runs doxygen and creates doxyerror
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile 2> ${CMAKE_BINARY_DIR}/DoxyError.full
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/..
		COMMENT "Generating API documentation with Doxygen")

	# filter the result of doxyerror
	if(UNIX)
		add_custom_command(
			TARGET doc
			POST_BUILD
			COMMAND grep -v "\"Detected potential recursive class relation\"" ${CMAKE_BINARY_DIR}/DoxyError.full 
			| grep -v "\"Member type (typedef)\"" 
			| grep -v "\"no uniquely matching class member\""
			| grep . >${CMAKE_BINARY_DIR}/DoxyError
			DEPENDS "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
			COMMENT "Filtering Doxygen errors and warnings"
		)
	endif(UNIX)

	# documentation installation
	if(NIHU_ENABLE_DOC_INSTALL)
		if(NOT (${NIHU_INSTALL_PATH} MATCHES ${CMAKE_BINARY_DIR}))
			message(STATUS "Documentation will be installed into ${CMAKE_INSTALL_PREFIX}/doc")
			install(DIRECTORY ${NIHU_HTML_DOC_DIR} DESTINATION doc)
		endif()
	endif(NIHU_ENABLE_DOC_INSTALL)

endif(DOXYGEN_FOUND)

### Installation section

# Select all hpp files for installation
set(NIHU_HPP_DIRECTORIES 
	"core" ; "interface" ; "library" ; "tmp" ; "util"; "aca")

foreach(HPP_DIRECTORY ${NIHU_HPP_DIRECTORIES})
	install(DIRECTORY ${HPP_DIRECTORY} DESTINATION include FILES_MATCHING PATTERN "*.hpp")
endforeach(HPP_DIRECTORY)

# Matlab installation section
set(NIHU_MATLAB_DIRECTORIES
	"analytic" ; "compatibility" ; "meshing")

foreach(MATLAB_DIRECTORY ${NIHU_MATLAB_DIRECTORIES})
	install(DIRECTORY "matlab/${MATLAB_DIRECTORY}" DESTINATION matlab FILES_MATCHING PATTERN "*.m")
endforeach(MATLAB_DIRECTORY)

# Installation rule for matlab demos
install(DIRECTORY "matlab/nihudemos" DESTINATION matlab)

# Create installation rule for matlab install script
install(FILES "matlab/install.m" DESTINATION matlab)

