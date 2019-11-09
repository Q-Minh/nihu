# NiHuLibraries.cmake 
# Library targets of NiHu

# Library suffices
if(WIN32)
	set(NIHU_STATIC_LIB_SUFFIX ".lib")
else()
	set(NIHU_STATIC_LIB_SUFFIX ".a")
endif()
	
### Library targets
# Create directory for the lib output
set(NIHU_LIB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/lib")
file(MAKE_DIRECTORY ${NIHU_LIB_OUTPUT_DIR})
# Add the library output directory to link 
link_directories(${NIHU_LIB_OUTPUT_DIR} ${Boost_LIBRARY_DIR_RELEASE})
if (DEFINED NIHU_BOOST_LIB_PATH)
	link_directories(${NIHU_BOOST_LIB_PATH})
endif()

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
if(TARGET "libfftw3-3")
	add_dependencies("nihu_static_lib" "libfftw3-3")
endif()

if(NOT DEFINED NIHU_DISABLE_FMM)
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
	if (TARGET "libfftw3-3")
		add_dependencies("nihu_fmm_static_lib" "libfftw3-3")
	endif()
endif()