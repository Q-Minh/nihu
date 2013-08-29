# Setup additional include directories
include_directories("${CMAKE_SOURCE_DIR}/util")

# Find all mex test sources
file(GLOB MEX_TEST_SOURCES *.mex.cpp)

# Create executable for all the unit tests
foreach (test_source ${MEX_TEST_SOURCES})
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${test_source})
	# Construct test name based on file name
	string(REPLACE ".mex.cpp" "" test_name ${local_source})
	# Construct the name of main m file
	set(test_mfile "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}_test.m")

	if(NOT NIHU_FORCE_MEX_COMPILER)
		# add the test as a shared library
		add_library(${test_name} SHARED ${local_source})
		# remove the "lib" prefix
		set_target_properties(${test_name} PROPERTIES 
			PREFIX "" 
			SUFFIX ${MATLAB_MEXEXT}
			COMPILE_FLAGS ${MEX_CXX_FLAGS}
			LINK_FLAGS ${MEX_SHARED_LINKER_FLAGS}
		)
	
	else(NOT NIHU_FORCE_MEX_COMPILER)

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
			COMMENT "Executing MEX for ${local_source}"
		)
	endif(NOT NIHU_FORCE_MEX_COMPILER)

	# copy the test m file
	add_custom_command(
		TARGET ${test_name} 
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy ${test_mfile} "./"
	)

	# create the script that executes matlab
	add_custom_command(
		TARGET ${test_name}
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -Dmat_root="${MATLAB_ROOT}" -Drun_file="run_${test_name}_test" -Dm_test_file="${test_name}_test.m" -P "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
		DEPENDS "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
	)

	if(NIHU_ENABLE_RUN_MATLAB_TESTS)
		# add the test
		if(WIN32)
			add_test(${test_name} "run_${test_name}_test.bat")
		elseif(WIN32)
			add_test(${test_name} "run_${test_name}_test")
		endif(WIN32)
	endif(NIHU_ENABLE_RUN_MATLAB_TESTS)

endforeach(test_source)
