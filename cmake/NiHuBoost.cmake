# NiHuBoost.cmake 
# Boost related

# Variables read by this script
#   NIHU_BOOST_PATH      Boost path specified by user 
#   NIHU_BOOST_LIB_PATH  Boost compiled libraries path (win only)
#
# Variables set by this script
#   NIHU_BOOST_LIBNAME_POSTFIX

# Process the boost path
message(STATUS "Looking for Boost ...")
if (DEFINED NIHU_BOOST_PATH)
    message(STATUS "\tBoost path (user input): " ${NIHU_BOOST_PATH})
	set(BOOST_ROOT ${NIHU_BOOST_PATH})
endif()

# Set variables for 
set(Boost_USE_STATIC_LIBS       OFF)  # find dynamic libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)  # use multi-threaded
set(Boost_USE_STATIC_RUNTIME    OFF)

# Find the boost package
# Boost >= 1.67 is required for advanced math functions, such as Chebyshev
if(WIN32)
	# On windows boost components are not found by FindBoost.cmake
	find_package(Boost 1.67 REQUIRED)
elseif(UNIX)
	find_package(Boost 1.67 REQUIRED COMPONENTS program_options)
endif()

if(BOOST_FOUND)
	message(STATUS "\tBoost was found, version: ${Boost_VERSION_STRING}")
	message(STATUS "\tBoost include directory: ${Boost_INCLUDE_DIRS}")
	if(NOT WIN32)
		message(STATUS "\tBoost components: ${Boost_LIBRARIES}")
		if(DEFINED Boost_PROGRAM_OPTIONS_LIBRARY)
			message(STATUS "\t Boost program options lib: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
		endif()
	endif()
else()
	message(STATUS "Boost was not found")
endif()

if(WIN32)
    if(MINGW)
		# No autolink for gcc 
		# Set postfix tag for linking libraries with gcc 
		set(NIHU_BOOST_LIBNAME_POSTFIX "mgw${GCC_MAJOR}${GCC_MINOR}-mt-x${NIHU_SYS_BITS}-${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}")
		message(STATUS "\tBoost library name postfix: -${NIHU_BOOST_LIBNAME_POSTFIX}")
	endif()
	# Check if library directory specified
	if(NOT DEFINED NIHU_BOOST_LIB_PATH)
	    message(WARNING "\tBoost library path was not set, link errors may occur")
		message(WARNING "\tTo avoid link errors specify NIHU_BOOST_LIB_PATH")
	endif()
endif()
