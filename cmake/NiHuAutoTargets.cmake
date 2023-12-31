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
	if("${target_mex_name}" STREQUAL "${target_name}")
		add_executable(${target_name} ${local_source})
		target_link_libraries(${target_name} ${NIHU_LINK_LIBRARIES} ${NIHU_FMM_LINK_LIBRARIES} ${FFTW3_LIBRARIES})
		if (WIN32 AND MINGW)
			target_link_libraries(${target_name}  "boost_program_options-${NIHU_BOOST_LIBNAME_POSTFIX}")
		elseif(DEFINED Boost_PROGRAM_OPTIONS_LIBRARY)
			target_link_libraries(${target_name} ${Boost_PROGRAM_OPTIONS_LIBRARY})
		endif()
		# Add installation
		install(TARGETS ${target_name} DESTINATION ${current_dir})
	elseif(NIHU_BUILD_MEX)
		# Do not use the mex compiler
		if(NOT NIHU_MATLAB_FORCE_MEX_COMPILER)
			# add the test as a shared library
			add_library(${target_mex_name} SHARED ${local_source})
			#target_link_libraries(${target_mex_name} ${NIHU_LINK_LIBRARIES_DYN}) 
			target_link_libraries(${target_mex_name} ${NIHU_LINK_LIBRARIES} ${NIHU_FMM_LINK_LIBRARIES} ${FFTW3_LIBRARIES} ${MEX_LINK_LIBRARIES})
			target_link_directories(${target_mex_name} PUBLIC ${MEX_LINK_DIRECTORIES})
			# remove the "lib" prefix
			set_target_properties(${target_mex_name} PROPERTIES 
				PREFIX "" 
				SUFFIX "${MATLAB_MEXEXT}"
				COMPILE_FLAGS "${MEX_CXX_FLAGS}"
				LINK_FLAGS "${MEX_SHARED_LINKER_FLAGS}"
			)
		# Use the mex compiler
		else()
			ADD_CUSTOM_TARGET (${target_mex_name} ALL)
			ADD_CUSTOM_COMMAND(
				TARGET    ${target_mex_name}
				COMMAND   ${MATLAB_MEX}
				ARGS      -O 
					CXXFLAGS=\""${CMAKE_CXX_FLAGS} -fPIC"\"
					-I"${MATLAB_ROOT}/extern/include" 
					-I"${EIGEN_INCLUDE_DIR}"
					-I"${CMAKE_SOURCE_DIR}" 
					"${CMAKE_CURRENT_SOURCE_DIR}/${local_source}"
					"${NIHU_COMMON_LIBRARIES}"
					-o "${test_name}"
				COMMENT "Executing MEX for ${local_source}"
			)
		endif()

		# Add installation
		install(TARGETS ${target_mex_name} DESTINATION ${current_dir})
	
		# Check if a corresponding .m file is present
		if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_mex_name}_app.m")
			set(target_mfile "${target_mex_name}_app.m")
		elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target_mex_name}_test.m")
			set(target_mfile "${target_mex_name}_test.m")
		else()
			unset(target_mfile)
		endif()

		# if there is a tester file
		if(target_mfile)
			# Copy the target m file to the build directory
			ADD_CUSTOM_COMMAND(
				TARGET ${target_mex_name}
				POST_BUILD
				COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${target_mfile}" "${CMAKE_CURRENT_BINARY_DIR}/"
				COMMENT "Copying .m file ${target_mfile} for target ${target_mex_name}${MATLAB_MEXEXT}"
			)
			# Copy the target .m file to the installation directory
			install(FILES ${target_mfile} DESTINATION ${current_dir})
		endif()
	endif()
endforeach()
