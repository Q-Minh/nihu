# NiHuMatlab.cmake
# Matlab related part of cmake for NiHu

#
# NIHU_MATLAB_FORCE_MEX_COMPILER
#   Force the usage of Matlab's mex compiler

message(STATUS "Looking for Matlab ...")

# look for environment variable MATLAB_ROOT
if(DEFINED NIHU_MATLAB_PATH)
	message(STATUS "\tMatlab path (user input): ${NIHU_MATLAB_PATH}")
	set(MATLAB_ROOT "${NIHU_MATLAB_PATH}")
elseif(DEFINED ENV{MATLAB_ROOT})
	message(STATUS "\tMatlab root (env. var.): $ENV{MATLAB_ROOT}")
	set(MATLAB_ROOT "$ENV{MATLAB_ROOT}")
endif()

# look for matlabmex pacakge
find_package (MatlabMex)

# check if matlab mex is found
if(MATLABMEX_FOUND)
 	message(STATUS "\tMatlab MEX found, enabling build of mex files")
 	set(NIHU_BUILD_MEX 1)
	# force usage of matlab's mex comipler
	if(NIHU_MATLAB_FORCE_MEX_COMPILER)
		message(STATUS "\tForcing Matlab's MEX compiler for mex files")
		message(STATUS "\t\tMEX compiler path: ${MATLAB_MEX}")
	else()
		message(STATUS "\tWill use the ${CMAKE_CXX_COMPILER_ID} compiler for building mex files")
		message(STATUS "\t\tCompiler path: ${CMAKE_CXX_COMPILER}")
	endif()
	# set the flags
	if(WIN32)
		set(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		set(MEX_LINK_DIRECTORIES 
			"${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}"
			"${MATLAB_ROOT}/extern/lib/win${NIHU_SYS_BITS}/microsoft"
			"${CMAKE_BINARY_DIR}/lib"
		)
		if(MSVC)
			# Correct for lib prefix in case of MSVC
			set(MEX_SHARED_LINKER_FLAGS "/export:mexFunction")
			set(MEX_LINK_LIBRARIES "libmex" "libmx" "libmat")
		else()
			set(MEX_SHARED_LINKER_FLAGS "-shared ${CMAKE_SHARED_LINKER_FLAGS} -m${NIHU_SYS_BITS}")
			set(MEX_LINK_LIBRARIES "mex" "mx" "mat")
		endif ()
	else()
		set(MEX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMATLAB_MEX_FILE")
		set(MEX_SHARED_LINKER_FLAGS "-L\"${CMAKE_BINARY_DIR}/lib\"")
		set(MEX_LINK_DIRECTORIES "${CMAKE_BINARY_DIR}/lib")
		set(MEX_LINK_LIBRARIES "")
	endif()
	# Set the extension for mex files
	set(MEX_SHARED_LIBRARY_SUFFIX "${MATLAB_MEXEXT}")
	# Add the include directory
	include_directories("${MATLAB_ROOT}/extern/include")
else()
	# if matlab mex was not found
	message(STATUS "\tMatlab MEX not found, disabling build of mex files")
	set(NIHU_BUILD_MEX 0)
endif()
