file(GLOB TEST_SOURCES *.cpp)

file(RELATIVE_PATH current_dir ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Create executable for all the unit tests
foreach (test_source ${TEST_SOURCES})
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${test_source})
	# Construct test name based on file name
	string(REPLACE ".cpp" "" test_name ${local_source})
	# Add the test executable
	add_executable(${test_name} ${local_source})
	# Add the test
	add_test(${test_name} ${test_name})
	# Install target
	if(NIHU_ENABLE_TEST_INSTALL)
		install(TARGETS ${test_name} DESTINATION ${current_dir})
	endif(NIHU_ENABLE_TEST_INSTALL)
endforeach(test_source)