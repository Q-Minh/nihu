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

### Linker configuration
if(WIN32)
	SET(CMAKE_EXE_LINKER_FLAGS "-MD ${CMAKE_EXE_LINKER_FLAGS} -s")
	SET(CMAKE_SHARED_LINKER_FLAGS "-MD ${CMAKE_SHARED_LINKER_FLAGS} -s")
endif(WIN32)

### EIGEN section

# Eigen as headers is default in windows mode
if(WIN32)
	# Disable eigen project install on windows
	set(NIHU_EIGEN_AS_PROJECT 0)
else(WIN32)
	# Turn of eigen as project as default on other operating systems
	if(NOT DEFINED NIHU_EIGEN_AS_PROJECT)
		set(NIHU_EIGEN_AS_PROJECT 0)
	endif(NOT DEFINED NIHU_EIGEN_AS_PROJECT)
endif(WIN32)

# Perform eigen installation
include(NiHuEigen)

### MATLAB Section

# look for environment variable MATLAB_ROOT
if (DEFINED ENV{MATLAB_ROOT})
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
elseif(DEFINED NIHU_MATLAB_ROOT)
	set(MATLAB_ROOT "${NIHU_MATLAB_ROOT}")
endif(DEFINED ENV{MATLAB_ROOT})

# look for matlabmex pacakge
find_package (MatlabMex)

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
	# set the flags
	if(WIN32)
		SET(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		SET(MEX_SHARED_LINKER_FLAGS "-shared ${CMAKE_SHARED_LINKER_FLAGS} -L\"${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}\" -L\"${MATLAB_ROOT}/extern/lib/win${NIHU_SYS_BITS}/microsoft\" -m${NIHU_SYS_BITS} -lmex -lmx -lmat")
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
	if(DEFINED NIHU_HTML_DOC_PATH)
	else(DEFINED NIHU_HTML_DOC_PATH)
		set(NIHU_HTML_DOC_PATH "${CMAKE_BINARY_DIR}/doc/html")
		file(MAKE_DIRECTORY ${NIHU_HTML_DOC_PATH})
	endif(DEFINED NIHU_HTML_DOC_PATH)
	message(STATUS "html documentation path: ${NIHU_HTML_DOC_PATH}")
	
	if(NOT NIHU_MATHJAX_DISABLE)
		message(STATUS "Enabling MathJax in Doxygen documentation")
		set(NIHU_MATHJAX_USE "YES")
		if(DEFINED NIHU_MATHJAX_PATH)
			file(RELATIVE_PATH NIHU_MATHJAX_RELPATH "${NIHU_HTML_DOC_PATH}" "${NIHU_MATHJAX_PATH}")
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
			install(DIRECTORY ${NIHU_HTML_DOC_PATH} DESTINATION doc)
		endif()
	endif(NIHU_ENABLE_DOC_INSTALL)

endif(DOXYGEN_FOUND)

### Installation section

# Select all hpp files for installation
set(NIHU_HPP_DIRECTORIES
	"bem"; "library"; "tmp"; "util")

foreach(HPP_DIRECTORY ${NIHU_HPP_DIRECTORIES})
	install(DIRECTORY ${HPP_DIRECTORY} DESTINATION include FILES_MATCHING PATTERN "*.hpp")
endforeach(HPP_DIRECTORY)