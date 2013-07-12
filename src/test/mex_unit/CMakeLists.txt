# Find all mex test sources
file(GLOB MEX_TEST_SOURCES *.mex.cpp)

# Create executable for all the unit tests
foreach (test_source ${MEX_TEST_SOURCES})
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${test_source})
	# Construct test name based on file name
	string(REPLACE ".mex.cpp" "" test_name ${local_source})

	ADD_CUSTOM_TARGET (${test_name} ALL)

	ADD_CUSTOM_COMMAND(
	TARGET    ${test_name}
	COMMAND   ${MATLAB_MEX}
	ARGS      -O 
		CXXFLAGS=\""${CMAKE_CXX_FLAGS} -fPIC"\"
		-I"${MATLAB_ROOT}/extern/include" 
		-I"${EIGEN_INCLUDE_DIR}"
		-I"${CMAKE_SOURCE_DIR}" 
		-I"${CMAKE_SOURCE_DIR}/util" 
		"${CMAKE_CURRENT_SOURCE_DIR}/${local_source}"
		-o "${test_name}"
	COMMENT   "Executing MEX for ${local_source}"
	)

	# TODO: add test runner
endforeach(test_source)