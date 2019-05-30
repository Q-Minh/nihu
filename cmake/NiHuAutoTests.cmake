# Find all mex test sources
file(GLOB MEX_TEST_SOURCES *.cpp)

# Remove unnecessary sources
list(REMOVE_ITEM MEX_TEST_SOURCES "${CMAKE_SOURCE_DIR}/test/core_unit/element_test.cpp")
list(REMOVE_ITEM MEX_TEST_SOURCES "${CMAKE_SOURCE_DIR}/test/core_unit/field_test.cpp")

file(RELATIVE_PATH current_dir ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Create executable for all the unit tests
foreach (test_source ${MEX_TEST_SOURCES})
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${test_source})
	# Construct test name based on file name
	string(REPLACE ".cpp" "" test_name ${local_source})
	string(REPLACE ".mex" "" test_mex_name ${test_name})
	
	# Cpp tests
	if(${test_name} STREQUAL ${test_mex_name})
		# Add the test executable
		if(${test_name} STREQUAL "core_test")
			add_executable(${test_name} ${local_source} 
				"element_test.cpp"
				"field_test.cpp")
		else()
			add_executable(${test_name} ${local_source})
		endif()
		target_link_libraries(${test_name} ${NIHU_LINK_LIBRARIES})
		# Add the test
		add_test(${test_name} ${test_name})
		# Install target
		if(NIHU_ENABLE_TEST_INSTALL)
			install(TARGETS ${test_name} DESTINATION ${current_dir})
		endif(NIHU_ENABLE_TEST_INSTALL)
	# MATLAB Tests
	elseif(NIHU_BUILD_MEX)
		# Find corresponding test file
		set(test_mfile "${CMAKE_CURRENT_SOURCE_DIR}/${test_mex_name}_test.m")

		if(NOT NIHU_MATLAB_FORCE_MEX_COMPILER)
			# add the test as a shared library
			add_library(${test_mex_name} SHARED ${local_source})
			#target_link_libraries(${test_mex_name} ${NIHU_LINK_LIBRARIES_DYN}) 
			target_link_libraries(${test_mex_name} ${NIHU_LINK_LIBRARIES} ${MEX_LINK_LIBRARIES}) 
			target_link_directories(${test_mex_name} PUBLIC ${MEX_LINK_DIRECTORIES})
			
			# remove the "lib" prefix
			set_target_properties(${test_mex_name} PROPERTIES 
				PREFIX "" 
				SUFFIX "${MATLAB_MEXEXT}"
				COMPILE_FLAGS "${MEX_CXX_FLAGS}"
				LINK_FLAGS "${MEX_SHARED_LINKER_FLAGS}"
			)
		
		else(NOT NIHU_MATLAB_FORCE_MEX_COMPILER)
			ADD_CUSTOM_TARGET (${test_mex_name} ALL)
			ADD_CUSTOM_COMMAND(
				TARGET    ${test_mex_name}
				COMMAND   ${MATLAB_MEX}
				ARGS      -O 
					CXXFLAGS=\""${CMAKE_CXX_FLAGS} -fPIC"\"
					-I"${MATLAB_ROOT}/extern/include" 
					-I"${EIGEN_INCLUDE_DIR}"
					-I"${CMAKE_SOURCE_DIR}" 
					-I"${CMAKE_SOURCE_DIR}/util" 
					"${CMAKE_CURRENT_SOURCE_DIR}/${local_source}"
					"${NIHU_COMMON_LIBRARIES}"
					-o "${test_name}"
				COMMENT "Executing MEX for ${local_source}"
			)
		endif(NOT NIHU_MATLAB_FORCE_MEX_COMPILER)

		# copy the test m file
		add_custom_command(
			TARGET ${test_mex_name} 
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${test_mfile} "${CMAKE_CURRENT_BINARY_DIR}/"
			COMMENT "Copying test .m file ${test_mfile} for target ${test_mex_name}${MATLAB_MEXEXT}"
		)

		# create the script that executes matlab
		add_custom_command(
			TARGET ${test_mex_name}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -Dmat_root="${MATLAB_ROOT}" -Drun_file="run_${test_name}_test" -Dm_test_file="${test_name}_test.m" -P "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
			DEPENDS "${CMAKE_MODULE_PATH}/WriteMatlabTestRunner.cmake"
			COMMENT "Creating run script for test .m file ${test_mfile}"
		)

		if(NIHU_ENABLE_RUN_MATLAB_TESTS)
			# add the test
			if(WIN32)
				add_test(${test_mex_name} "run_${test_name}_test.bat")
			elseif(WIN32)
				add_test(${test_mex_name} "run_${test_name}_test")
			endif(WIN32)
		endif(NIHU_ENABLE_RUN_MATLAB_TESTS)

		if(NIHU_ENABLE_TEST_INSTALL)
			install(TARGETS ${test_mex_name} DESTINATION ${current_dir})
			install(FILES ${test_mfile} DESTINATION ${current_dir})
		endif(NIHU_ENABLE_TEST_INSTALL)
	endif()

endforeach(test_source)
