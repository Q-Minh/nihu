file(GLOB CPP_SOURCES *.cpp)

file(RELATIVE_PATH current_dir ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Create executable for all the unit tests
foreach (cpp_source ${CPP_SOURCES})
	# get current directory
	# get_filename_component(local_source ${test_source} RELATIVE)
	file(RELATIVE_PATH local_source ${CMAKE_CURRENT_SOURCE_DIR} ${cpp_source})
	# Construct test name based on file name
	string(REPLACE ".cpp" "" target_name ${local_source})
	string(REPLACE ".mex" "" target_mex_name ${target_name})

	# Add the executable / or mex file
	if(${target_mex_name} MATCHES ${target_name})
		add_executable(${target_name} ${local_source})
		# Add installation
		install(TARGETS ${target_name} DESTINATION ${current_dir})
	else()
		# Do not use the mex compiler
		if(NOT NIHU_FORCE_MEX_COMPILER)
			# add the test as a shared library
			add_library(${target_mex_name} SHARED ${local_source})
			# remove the "lib" prefix
			set_target_properties(${target_mex_name} PROPERTIES 
				PREFIX "" 
				SUFFIX "${MATLAB_MEXEXT}"
				COMPILE_FLAGS "${MEX_CXX_FLAGS}"
				LINK_FLAGS "${MEX_SHARED_LINKER_FLAGS}"
			)
		# Use the mex compiler
		else(NOT NIHU_FORCE_MEX_COMPILER)
			ADD_CUSTOM_TARGET (${target_mex_name} ALL)
			ADD_CUSTOM_COMMAND(
				TARGET    ${target_mex_name}
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

		# Add installation
		install(TARGETS ${target_mex_name} DESTINATION ${current_dir})
	endif()

	
	
endforeach(cpp_source)