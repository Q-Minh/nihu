# cmake file for the NiHu project

###
### cmake and project section
###

# Set the minimum required cmake version
cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)

# Set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# Set the project name
project (NiHu)

set (NIHU_VERSION_MAJOR 2)
set (NIHU_VERSION_MINOR 0)

# Find the number of system bits
include(Sysbits)
get_system_bits(NIHU_SYS_BITS)
message(STATUS "${NIHU_SYS_BITS}-bit system detected")

# Set c++11 as requirement
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

### Linker configuration
if(WIN32)
	SET(CMAKE_EXE_LINKER_FLAGS "-MD ${CMAKE_EXE_LINKER_FLAGS} -s")
	SET(CMAKE_SHARED_LINKER_FLAGS "-MD ${CMAKE_SHARED_LINKER_FLAGS} -s")
endif(WIN32)

###
### ThirdParty 
###

# Set the 3rdparty dir
set(NIHU_THIRDPARTY_DIR "${CMAKE_SOURCE_DIR}/ThirdParty")

### EIGEN section
include(NiHuEigen)

### Boost section
# Process the boost path
if (DEFINED NIHU_BOOST_PATH)
	set(BOOST_ROOT ${NIHU_BOOST_PATH})
endif ()

# Find the boost package
find_package(Boost 1.67 REQUIRED)
if (BOOST_FOUND)
	message(STATUS "Boost was found, include directory: ${Boost_INCLUDE_DIR}")
else ()
	message(STATUS "Boost was not found")
endif ()

### FFTW section 
# Process fftw path
if (DEFINED NIHU_FFTW_PATH)
	set(FFTW3_DIRECTORY ${NIHU_FFTW_PATH})
endif ()

find_package(FFTW REQUIRED)
if (FFTW3_FOUND)
	message(STATUS "fftw3 was found, include directory: ${FFTW3_INCLUDE_DIRS}")
endif ()

# Add FFTW include directory
include_directories(${FFTW3_INCLUDE_DIRS})

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
		set(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		set(MEX_SHARED_LINKER_FLAGS "-L\"${CMAKE_BINARY_DIR}/lib\"")
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



# Setup compiler flags for gcc
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pedantic -O3 -Wall \
		-Wno-deprecated-declarations")
	set(NIHU_FPIC_FLAG "-fPIC")
	# Find out gcc version
	execute_process(
		COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
		OUTPUT_VARIABLE GCC_VERSION
	)
    # Add -Wno-misleading-indentation for gcc >= 6
    # Note: this disables Eigen-generated warnings
    if (GCC_VERSION VERSION_GREATER 6.0 OR GCC_VERSION VERSION_EQUAL 6.0)
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} \
			-Wno-misleading-indentation")
    endif (GCC_VERSION VERSION_GREATER 6.0 OR GCC_VERSION VERSION_EQUAL 6.0)
    # Add 
    # Note: this disables Eigen-generated warnings
    if (GCC_VERSION VERSION_GREATER 7.0 OR GCC_VERSION VERSION_EQUAL 7.0)
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} \
			-Wno-int-in-bool-context \
			-Wno-maybe-uninitialized \
			-Wno-uninitialized")
    endif (GCC_VERSION VERSION_GREATER 7.0 OR GCC_VERSION VERSION_EQUAL 7.0)
endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

### Includes
# setup global include directories
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${EIGEN_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIR}")

### Libraries
if(WIN32)
	set(NIHU_STATIC_LIB_SUFFIX ".lib")
else(WIN32)
	set(NIHU_STATIC_LIB_SUFFIX ".a")
endif(WIN32)
	
### Library targets
# Create directory for the lib output
set(NIHU_LIB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/lib")
file(MAKE_DIRECTORY ${NIHU_LIB_OUTPUT_DIR})
# Add the library output directory to link 
link_directories(${NIHU_LIB_OUTPUT_DIR})

# Look for common library source files
FILE(GLOB NIHU_COMMON_LIBRARIES 
	"${CMAKE_SOURCE_DIR}/library/lib_*.cpp")

# Add the static NiHu library
add_library("nihu_static_lib" STATIC ${NIHU_COMMON_LIBRARIES})
set_target_properties("nihu_static_lib"
	PROPERTIES
	OUTPUT_NAME "nihu"
	COMPILE_FLAGS "${NIHU_FPIC_FLAG}"
	PREFIX "" 
	SUFFIX ${NIHU_STATIC_LIB_SUFFIX}
	ARCHIVE_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
	LIBRARY_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
)
set(NIHU_LINK_LIBRARIES ${NIHU_LINK_LIBRARIES}; "nihu_static_lib")

if (0)
# Add the dynamic NiHu library
add_library("nihu_shared_lib" SHARED ${NIHU_COMMON_LIBRARIES})
set_target_properties("nihu_shared_lib"
	PROPERTIES 
	COMPILE_FLAGS "${NIHU_FPIC_FLAG}"
	OUTPUT_NAME "nihu"
	PREFIX ""
	ARCHIVE_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
	LIBRARY_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
)
set(NIHU_LINK_LIBRARIES_DYN ${NIHU_LINK_LIBRARIES_DYN}; "nihu_shared_lib") 
endif (0)

# Add the static FMM library 
# Look for common library source files
FILE(GLOB NIHU_FMM_LIBRARY_SOURCES
	"${CMAKE_SOURCE_DIR}/fmm/*.cpp")
add_library("nihu_fmm_static_lib" STATIC ${NIHU_FMM_LIBRARY_SOURCES})
set_target_properties("nihu_fmm_static_lib"
	PROPERTIES
	OUTPUT_NAME "nihu_fmm"
	COMPILE_FLAGS "${NIHU_FPIC_FLAG}"
	PREFIX "" 
	SUFFIX ${NIHU_STATIC_LIB_SUFFIX}
	ARCHIVE_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
	LIBRARY_OUTPUT_DIRECTORY ${NIHU_LIB_OUTPUT_DIR}
)
set(NIHU_FMM_LINK_LIBRARIES ${NIHU_FMM_LINK_LIBRARIES}; "nihu_fmm_static_lib")


### Subdirectories
# Tutorials directory
if(NOT NIHU_DISABLE_TUTORIALS)
	add_subdirectory(tutorial)
endif()

### Test section
if(NIHU_ENABLE_TESTING)
	# Check for gtest module
	include (NiHuGTest)
	# Check if gtest directories are found
	if(NOT DEFINED GTEST_INCLUDE_DIRS)
		message(STATUS "GTest not found, or could not be installed.")
		message(STATUS "Build of gtest tests are disabled.")
	else ()
		message(STATUS "GTest headers will be imported from: ${GTEST_INCLUDE_DIRS}")
		include_directories("${GTEST_INCLUDE_DIRS}")
		include_directories("${GTEST_MAIN_DIR}")
		add_subdirectory(gtest)
	endif()
	
	# Check if install 
	if(NIHU_ENABLE_TEST_INSTALL)
		if(${NIHU_INSTALL_DIR} MATCHES ${CMAKE_BINARY_DIR})
			set(NIHU_ENABLE_TEST_INSTALL 0)
			message(STATUS "Installation of tests disabled (build dir is same as install dir)")
		endif()
	endif(NIHU_ENABLE_TEST_INSTALL)
	# Enable testing globally
	enable_testing()
	# Check if automatic build of tests are disabled
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
### Setup directories for installation
if(NOT DEFINED NIHU_INSTALL_DIR)
	set(NIHU_INSTALL_DIR ${CMAKE_BINARY_DIR})
endif(NOT DEFINED NIHU_INSTALL_DIR)

set(CMAKE_INSTALL_PREFIX ${NIHU_INSTALL_DIR})

# Select all hpp files for installation
set(NIHU_HPP_DIRECTORIES 
	"core" ; "interface" ; "library" ; "tmp" ; "util" ; "aca")

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

