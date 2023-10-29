# CMakeLists.txt
# cmake file for the NiHu project
# 
# This file and its submodules follow the coding syntax defined here:
#   https://community.kde.org/Policies/CMake_Coding_Style

# NOTE: There is a problem of Win10 + Matlab R2018a / symbolic math toolbox

### MATLAB section
if (NIHU_ENABLE_MATLAB)
	include(${CMAKE_CURRENT_LIST_DIR}/../cmake/NiHuMatlab.cmake)
endif()

### GTest section
if (NIHU_ENABLE_TESTING)
	# Check for gtest module
	include(${CMAKE_CURRENT_LIST_DIR}/../cmake/NiHuGTest.cmake)
endif()
	
###
### Modules configuration 
###

### FMM module configuration
if (NIHU_ENABLE_FMM)
	include(${CMAKE_CURRENT_LIST_DIR}/../cmake/NiHuFMM.cmake)
endif()

###
### Library targets
###

add_library(nihu STATIC)
add_library(nihu::nihu ALIAS nihu)

set(_nihu_source_file_patterns "")
list(APPEND _nihu_source_file_patterns 
	"aca/*.cpp" 
	"core/*.cpp" 
	"cqm/*.cpp" 
	"library/*.cpp")
set(_nihu_header_file_patterns "")
list(APPEND _nihu_header_file_patterns 
	"aca/*.h" 
	"aca/*.hpp" 
	"core/*.h" 
	"core/*.hpp" 
	"cqm/*.h" 
	"cqm/*.hpp" 
	"library/*.h" 
	"library/*.hpp" 
	"tmp/*.h" 
	"tmp/*.hpp" 
	"util/*.h"
	"util/*.hpp"
	)

if (NIHU_ENABLE_FMM)
	list(APPEND _nihu_source_file_patterns "fmm/*.cpp")
	list(APPEND _nihu_header_file_patterns "fmm/*.h" "fmm/*.hpp")
endif()

file(
    GLOB_RECURSE _nihu_sources 
    LIST_DIRECTORIES OFF 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    CONFIGURE_DEPENDS 
    ${_nihu_source_file_patterns}
    )
file(
    GLOB_RECURSE _nihu_headers 
    LIST_DIRECTORIES OFF 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    CONFIGURE_DEPENDS 
    ${_nihu_header_file_patterns}
	)

target_sources(nihu 
PRIVATE
	${_nihu_sources}
	${_nihu_headers}
)

source_group(
TREE 
    "${CMAKE_CURRENT_SOURCE_DIR}" 
FILES 
    ${_nihu_sources} 
    ${_nihu_headers}
)

set_target_properties(nihu 
PROPERTIES 
    FOLDER nihu
    WINDOWS_EXPORT_ALL_SYMBOLS ON 
    DEBUG_POSTFIX "_d" 
    POSITION_INDEPENDENT_CODE ON
)

target_compile_definitions(nihu PUBLIC cxx_std_17)
target_compile_options(nihu PRIVATE /bigobj)
include(GNUInstallDirs)
target_include_directories(nihu 
PUBLIC 
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
)
target_link_libraries(nihu 
PUBLIC 
	Eigen3::Eigen 
	Boost::boost 
	FFTW3::fftw3
	OpenMP::OpenMP_CXX)

# Tutorials directory
if (NIHU_ENABLE_TUTORIALS)
	add_subdirectory(tutorial)
endif()

# Applications/examples
if (NIHU_ENABLE_APPLICATIONS)
	add_subdirectory(applications/helmholtz)
	add_subdirectory(applications/laplace)
	add_subdirectory(applications/stochastic)
endif()

### Test section
if (NIHU_ENABLE_TESTING)
	include(${CMAKE_CURRENT_LIST_DIR}/../cmake/NiHuTesting.cmake)
endif()

### Documentation section
if (NIHU_ENABLE_DOC)
	include(${CMAKE_CURRENT_LIST_DIR}/../cmake/NiHuDoc.cmake)
endif()
