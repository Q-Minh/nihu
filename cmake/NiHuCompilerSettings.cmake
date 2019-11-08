# Setup compiler flags for gcc
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pedantic -O3 -Wall -fopenmp \
		-Wno-deprecated-declarations")
	set(NIHU_FPIC_FLAG "-fPIC")
	# Find out gcc version
	execute_process(
		COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
		OUTPUT_VARIABLE GCC_VERSION
	)
	# Store the GCC_MAJOR and GCC_MINOR versions for later use
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
	list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
	list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
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
			-Wno-parentheses \
			-Wno-sign-compare \
			-Wno-uninitialized")
    endif (GCC_VERSION VERSION_GREATER 7.0 OR GCC_VERSION VERSION_EQUAL 7.0)
endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

# Setup compiler flags for MSVC
if (MSVC)
	MESSAGE(STATUS "Adding /bigobj compiler flag")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /openmp")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
endif (MSVC)
