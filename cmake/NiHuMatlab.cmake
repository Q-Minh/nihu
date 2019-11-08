# look for environment variable MATLAB_ROOT
if (DEFINED ENV{MATLAB_ROOT})
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
elseif(DEFINED NIHU_MATLAB_PATH)
	set(MATLAB_ROOT "${NIHU_MATLAB_PATH}")
endif(DEFINED ENV{MATLAB_ROOT})

# look for matlabmex pacakge
find_package (MatlabMex)

# check if matlab mex is found
if(MATLABMEX_FOUND)
 	message(STATUS "Matlab MEX found, enabling build of mex files")
 	set(NIHU_BUILD_MEX 1)
	# force usage of matlab's mex comipler
	if(NIHU_MATLAB_FORCE_MEX_COMPILER)
		message(STATUS "Forcing Matlab's MEX compiler (${MATLAB_MEX}) for mex files")
	else(NIHU_MATLAB_FORCE_MEX_COMPILER)
		message(STATUS "Will use the ${CMAKE_CXX_COMPILER_ID} compiler (${CMAKE_CXX_COMPILER}) for mex files")
	endif(NIHU_MATLAB_FORCE_MEX_COMPILER)
	# set the flags
	if(WIN32)
		SET(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		#SET(MEX_SHARED_LINKER_FLAGS "-shared ${CMAKE_SHARED_LINKER_FLAGS} -m${NIHU_SYS_BITS} -L\"${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}\" -L\"${MATLAB_ROOT}/extern/lib/win${NIHU_SYS_BITS}/microsoft\" -L\"${CMAKE_BINARY_DIR}/lib\"  -lmex -lmx -lmat")
		set(MEX_LINK_DIRECTORIES 
			"${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}"
			"${MATLAB_ROOT}/extern/lib/win${NIHU_SYS_BITS}/microsoft"
			"${CMAKE_BINARY_DIR}/lib"
		)
		if (MSVC)
			# Correct for lib prefix in case of MSVC
			set(MEX_SHARED_LINKER_FLAGS "/export:mexFunction")
			set(MEX_LINK_LIBRARIES "libmex" "libmx" "libmat")
		else ()
			set(MEX_SHARED_LINKER_FLAGS "-shared ${CMAKE_SHARED_LINKER_FLAGS} -m${NIHU_SYS_BITS}")
			set(MEX_LINK_LIBRARIES "mex" "mx" "mat")
		endif ()
	else(WIN32)
		set(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		set(MEX_SHARED_LINKER_FLAGS "-L\"${CMAKE_BINARY_DIR}/lib\"")
		set(MEX_LINK_DIRECTORIES "${CMAKE_BINARY_DIR}/lib")
		set(MEX_LINK_LIBRARIES "")
	endif(WIN32)
	# Set the extension for mex files
	set(MEX_SHARED_LIBRARY_SUFFIX "${MATLAB_MEXEXT}")
	# Add the include directory
	include_directories("${MATLAB_ROOT}/extern/include")
# if matlab mex is not found
else(MATLABMEX_FOUND)
	message(STATUS "Matlab MEX not found, disabling build of mex files")
	set(NIHU_BUILD_MEX 0)
endif(MATLABMEX_FOUND)