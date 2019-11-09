# This module looks for mex, the MATLAB compiler.
# The following variables are defined when the script completes:
#   MATLAB_MEX: location of mex compiler
#   MATLAB_ROOT: root of MATLAB installation
#   MATLABMEX_FOUND: 0 if not found, 1 if found
# 	MATLAB_MEXEXT: matlab mex extension

set(MATLAB_MEX "")
set(MATLABMEX_FOUND 0)

if(WIN32)
	set(MATLAB_MEXEXT .mexw)
	# Modifications by Ramon Casero and Tom Doel for Gerardus.
	# Updated by Peter Rucz, 2013.
	# Further updated by Peter Rucz, 2019.

	# Search for a version of Matlab available, starting from the most modern one
	# to older versions.
	if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
		foreach(MATVER 
			"9.6" "9.5" "9.4" "9.3" "9.2" "9.1" "9.0" 
			"8.6" "8.5" "8.4" "8.3" "8.2" "8.1" "8.0" 
			"7.14" "7.13" "7.12" "7.11" "7.10" "7.9" "7.8" "7.7" "7.6" "7.5" "7.4" 
		)
			if((NOT DEFINED MATLAB_ROOT)
				OR ("${MATLAB_ROOT}" STREQUAL "")
				OR ("${MATLAB_ROOT}" STREQUAL "/registry"))
				get_filename_component(MATLAB_ROOT
					"[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${MATVER};MATLABROOT]"
					ABSOLUTE)
				set(MATLAB_VERSION ${MATVER})
			endif()
		endforeach()
	endif()

	if(EXISTS "${MATLAB_ROOT}/bin/mex.bat")
		set(MATLAB_MEX "${MATLAB_ROOT}/bin/mex.bat")
	endif(EXISTS "${MATLAB_ROOT}/bin/mex.bat")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Check if this is a Mac.
	set(MATLAB_MEXEXT .maci)
	
	# This code is untested but taken from the older FindMatlab.cmake script as
	# well as the modifications by Ramon Casero and Tom Doel for Gerardus.
	set(LIBRARY_EXTENSION .dylib)

	# If this is a Mac and the attempts to find MATLAB_ROOT have so far failed,~
	# we look in the applications folder
	if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))

		# Search for a version of Matlab available, starting from the most modern
		# one to older versions
		foreach(MATVER 
			"R2019b" "R2019a" 
			"R2018b" "R2018a" "R2017b" "R2017a" "R2016b" "R2016a" 
			"R2015b" "R2015a" "R2014b" "R2014a" "R2013b" "R2013a" 
			"R2012b" "R2012a" "R2011b" "R2011a" "R2010b" "R2010a" 
			"R2009b" "R2009a" "R2008b"
		)
			if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
				if(EXISTS /Applications/MATLAB_${MATVER}.app)
					set(MATLAB_ROOT /Applications/MATLAB_${MATVER}.app)
				endif()
			endif()
		endforeach()
	endif()

	find_program(MATLAB_MEX
		mex
		PATHS
		${MATLAB_ROOT}/bin
		${MATLAB_ROOT}/bin/win${NIHU_SYS_BITS}
	)
elseif(UNIX)
	# On a Linux system.  The goal is to find MATLAB_ROOT.
	set(MATLAB_MEXEXT .mexa)
	set(LIBRARY_EXTENSION .so)

	if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
		execute_process(COMMAND "/usr/bin/locate" "-e" "--regex" "[Mm][Aa][Tt][Ll][Aa][Bb].*/bin/mex$"
			RESULT_VARIABLE locate_matlab_result
			OUTPUT_VARIABLE locate_matlab_output
		)

		# If locating is successful and paths are found
		if(NOT ${locate_matlab_result} AND NOT "${locate_matlab_output}" STREQUAL "")
			# Convert result into list
			string(REPLACE "\n" ";" locate_matlab_list ${locate_matlab_output})
			# Revert the order of the list
			list(REVERSE locate_matlab_list)
			list(GET locate_matlab_list 1 locate_matlab_path)
		
			# Get the path of matlab mex executable
			get_filename_component(MATLAB_MEX "${locate_matlab_path}" REALPATH)
			get_filename_component(MATLAB_BIN_ROOT "${MATLAB_MEX}" PATH)
			# Strip ./bin/.
			GET_FILENAME_COMPONENT(MATLAB_ROOT "${MATLAB_BIN_ROOT}" PATH)
		
		# If locating was not successful
		else(NOT ${locate_matlab_result} AND NOT "${locate_matlab_output}" STREQUAL "")
			find_program(MATLAB_MEX_POSSIBLE_LINK
				NAMES
				mex
				PATHS
				${MATLAB_ROOT}/bin
				/opt/matlab/*/bin
				/usr/local/matlab/*/bin
				$ENV{HOME}/matlab/*/bin
				NO_DEFAULT_PATH
			)
			# Check if the mex link is found
			if(NOT "${MATLAB_MEX_POSSIBLE_LINK}" STREQUAL "MATLAB_MEX_POSSIBLE_LINK-NOTFOUND" AND NOT "${MATLAB_MEX_POSSIBLE_LINK}" STREQUAL "")
				get_filename_component(MATLAB_MEX "${MATLAB_MEX_POSSIBLE_LINK}" REALPATH)
				get_filename_component(MATLAB_BIN_ROOT "${MATLAB_MEX}" PATH)
				# Strip ./bin/.
				get_filename_component(MATLAB_ROOT "${MATLAB_BIN_ROOT}" PATH)
			else()
				set(MATLAB_ROOT "")
				set(MATLAB_MEX "")
			endif()
		endif()
	else()
		# MATLAB_ROOT defined and not empty
		if(EXISTS "${MATLAB_ROOT}/bin/mex")
			set(MATLAB_MEX "${MATLAB_ROOT}/bin/mex")
		else()
			set(MATLAB_ROOT "")
			set(MATLAB_MEX "")
		endif()
	endif()
endif()

if("${MATLAB_MEX}" STREQUAL "")
	# Matlab mex not found
	message(STATUS "\tCould not find MATLAB mex compiler")
	message(STATUS "\t\tTry specifying NIHU_MATLAB_PATH if Matlab is installed")
else()
	set(MATLAB_MEXEXT "${MATLAB_MEXEXT}${NIHU_SYS_BITS}")
	message(STATUS "\tFound MATLAB mex compiler: ${MATLAB_MEX}")
	message(STATUS "\tMATLAB root: ${MATLAB_ROOT}")
	message(STATUS "\tMEX file extension is: ${MATLAB_MEXEXT}")
	set(MATLABMEX_FOUND 1)
endif()
