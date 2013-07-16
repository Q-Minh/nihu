add_subdirectory(tmp_unit)
add_subdirectory(core_unit)
add_subdirectory(util_unit)
add_subdirectory(library_unit)

if(NIHU_BUILD_MEX)
	add_subdirectory(mex_unit)
	add_subdirectory(matlab)
endif(NIHU_BUILD_MEX)

