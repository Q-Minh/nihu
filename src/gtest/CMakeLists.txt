# Find all mex test sources
file(GLOB TEST_SOURCES *.cpp)

set(gtest_main_name "test_gtest")

add_executable(${gtest_main_name} ${TEST_SOURCES} ${GTEST_MAIN_DIR}/src/gtest_main.cc ${GTEST_MAIN_DIR}/src/gtest-all.cc)

set_target_properties(${gtest_main_name}
	PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${GTEST_INCLUDE_DIRS}")
target_link_libraries(${gtest_main_name} ${NIHU_LINK_LIBRARIES} pthread) 