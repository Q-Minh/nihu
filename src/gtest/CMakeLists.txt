# TODO:
# 	Found error on Win7 + MinGW64 8.1: the built test_gtest.exe tries to
#	load a different version of libgomp.dll and libstdc++-6.dll copying
# 	the files from MinGW/bin beside the test_gtest.exe solves the problem.
#	This may be a path issue.
#	The path conflict seems to be coming from gstreamer install.

# Find all mex test sources
file(GLOB TEST_SOURCES *.cpp)

set(gtest_main_name "test_gtest")

add_executable(${gtest_main_name} ${TEST_SOURCES} ${GTEST_MAIN_DIR}/src/gtest-all.cc)
target_include_directories(${gtest_main_name} PUBLIC "${GTEST_INCLUDE_DIRS}" "${GTEST_MAIN_DIR}")

if(MSVC)
	target_link_libraries(${gtest_main_name} ${NIHU_LINK_LIBRARIES} ${NIHU_FMM_LINK_LIBRARIES} ${FFTW3_LIBRARIES})
else()
	if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
		# Add the -Wno-unused-variable flag for gcc to build without warnings
		set_target_properties(${gtest_main_name}
			PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
		target_link_libraries(${gtest_main_name} ${NIHU_LINK_LIBRARIES} ${NIHU_FMM_LINK_LIBRARIES} ${FFTW3_LIBRARIES} pthread)
	endif()
endif()
