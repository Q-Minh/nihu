file(GLOB CPP_SOURCES *.cpp)

file(RELATIVE_PATH current_dir ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Create executable for all the unit tests
foreach (cpp_source ${CPP_SOURCES})
	# get current directory
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${cpp_source})
	# Construct test name based on file name
	string(REPLACE ".cpp" "" test_name ${local_source})
	# Add the test executable
	add_executable(${test_name} ${local_source})
	# Add installation
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${test_name}" DESTINATION "bin/${current_dir}")
endforeach(cpp_source)