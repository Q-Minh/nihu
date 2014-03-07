# Test directories with cpp tests
set(cpp_test_dirs
	"core_unit"; "library_unit"; "stallone_unit"; "tmp_unit"; "util_unit")

# Add each test directory
foreach(cpp_test_dir ${cpp_test_dirs})
	# Check if directory exists
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}")
		if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}/CMakeLists.txt")
			add_subdirectory(${cpp_test_dir})
		else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}/CMakeLists.txt")
			message(WARNING "Could not add C++ test directory: ${cpp_test_dir}, no CMakeLists.txt present")
		endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}/CMakeLists.txt")
	# If no CMakeLists file
	else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}")
		message(WARNING "Could not add C++ test directory: ${cpp_test_dir}, directory not found")
	endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cpp_test_dir}")
endforeach(cpp_test_dir)

# Test directories with matlab tests 
set(matlab_test_dirs
	"mex_unit"; "matlab")

if(NIHU_BUILD_MEX)
	foreach(matlab_test_dir ${matlab_test_dirs})
		# Check if directory exists
		if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}")
			if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}/CMakeLists.txt")
				add_subdirectory(${matlab_test_dir})
			else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}/CMakeLists.txt")
				message(WARNING "Could not add Matlab test directory: ${cpp_test_dir}, no CMakeLists.txt present")
			endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}/CMakeLists.txt")
		# If no CMakeLists file
		else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}")
			message(WARNING "Could not add Matlab test directory: ${matlab_test_dir}, directory not found")
		endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${matlab_test_dir}")
	endforeach(matlab_test_dir)
endif(NIHU_BUILD_MEX)
