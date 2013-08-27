# Test directories with cpp tests
set(cpp_test_dirs
	"core_unit"; "library_unit"; "tmp_unit"; "util_unit"; "integral")

foreach(cpp_test_dir ${cpp_test_dirs})
	add_subdirectory(${cpp_test_dir})
endforeach(cpp_test_dir)

# Test directories with matlab tests 
set(matlab_test_dirs
	"mex_unit"; "matlab")

if(NIHU_BUILD_MEX)
	foreach(matlab_test_dir ${matlab_test_dirs})
		add_subdirectory(${matlab_test_dir})
	endforeach(matlab_test_dir)
endif(NIHU_BUILD_MEX)
