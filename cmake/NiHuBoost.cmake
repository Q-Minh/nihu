# Process the boost path
message(STATUS "Looking for Boost ...")
if (DEFINED NIHU_BOOST_PATH)
    message(STATUS "\tBoost path (user input): " ${NIHU_BOOST_PATH})
	set(BOOST_ROOT ${NIHU_BOOST_PATH})
endif ()

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)  # use multi-threaded
set(Boost_USE_STATIC_RUNTIME    OFF)

# Find the boost package
if (WIN32)
	find_package(Boost 1.67 REQUIRED)
else ()
	find_package(Boost 1.67 REQUIRED COMPONENTS math program_options)
endif ()
if (BOOST_FOUND)
	message(STATUS "\tBoost was found, version: ${Boost_VERSION_STRING}")
	message(STATUS "\tBoost include directory: ${Boost_INCLUDE_DIRS}")
	if (NOT WIN32)
		message(STATUS "\tBoost components: ${Boost_LIBRARIES}")
	endif()
else ()
	message(STATUS "Boost was not found")
endif ()

# No autolink for gcc set postfix
if (MINGW)
	# Set postfix tag for linking libraries with gcc 
	SET(NIHU_BOOST_LIBNAME_POSTFIX "mgw${GCC_MAJOR}${GCC_MINOR}-mt-x${NIHU_SYS_BITS}-${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}")
	message(STATUS "\tBoost library name postfix: -${NIHU_BOOST_LIBNAME_POSTFIX}")
endif()