# This module looks for mex, the MATLAB compiler.
# The following variables are defined when the script completes:
#   MATLAB_MEX: location of mex compiler
#   MATLAB_ROOT: root of MATLAB installation
#   MATLABMEX_FOUND: 0 if not found, 1 if found
# 	MATLAB_MEXEXT: matlab mex extension

SET(MATLAB_MEX "")
SET(MATLABMEX_FOUND 0)

IF(WIN32)
	SET(MATLAB_MEXEXT .mexw)
	# Modifications by Ramon Casero and Tom Doel for Gerardus.
	# Updated by Peter Rucz, 2013.
	# Further updated by Peter Rucz, 2019.

	# Search for a version of Matlab available, starting from the most modern one
	# to older versions.
	if((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
	FOREACH(MATVER 
		"9.6" "9.5" "9.4" "9.3" "9.2" "9.1" "9.0" 
		"8.6" "8.5" "8.4" "8.3" "8.2" "8.1" "8.0" 
		"7.14" "7.13" "7.12" "7.11" "7.10" "7.9" "7.8" "7.7" "7.6" "7.5" "7.4" 
	)
		IF((NOT DEFINED MATLAB_ROOT)
			OR ("${MATLAB_ROOT}" STREQUAL "")
			OR ("${MATLAB_ROOT}" STREQUAL "/registry"))
			GET_FILENAME_COMPONENT(MATLAB_ROOT
				"[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${MATVER};MATLABROOT]"
				ABSOLUTE)
			SET(MATLAB_VERSION ${MATVER})
		ENDIF((NOT DEFINED MATLAB_ROOT)
			OR ("${MATLAB_ROOT}" STREQUAL "")
			OR ("${MATLAB_ROOT}" STREQUAL "/registry"))
	ENDFOREACH(MATVER)
	endif()

	if(EXISTS "${MATLAB_ROOT}/bin/mex.bat")
		set(MATLAB_MEX "${MATLAB_ROOT}/bin/mex.bat")
	endif(EXISTS "${MATLAB_ROOT}/bin/mex.bat")

ELSE(WIN32 OR MSYS)
	SET(MATLAB_MEXEXT .maci)
	# Check if this is a Mac.
	IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		# This code is untested but taken from the older FindMatlab.cmake script as
		# well as the modifications by Ramon Casero and Tom Doel for Gerardus.

		SET(LIBRARY_EXTENSION .dylib)

		# If this is a Mac and the attempts to find MATLAB_ROOT have so far failed,~
		# we look in the applications folder
		IF((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))

			# Search for a version of Matlab available, starting from the most modern
			# one to older versions
			FOREACH(MATVER 
				"R2019a" 
				"R2018b" "R2018a" "R2017b" "R2017a" "R2016b" "R2016a" 
				"R2015b" "R2015a" "R2014b" "R2014a" "R2013b" "R2013a" 
				"R2012b" "R2012a" "R2011b" "R2011a" "R2010b" "R2010a" 
				"R2009b" "R2009a" "R2008b"
			)
				IF((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
					IF(EXISTS /Applications/MATLAB_${MATVER}.app)
						SET(MATLAB_ROOT /Applications/MATLAB_${MATVER}.app)
					ENDIF(EXISTS /Applications/MATLAB_${MATVER}.app)
				ENDIF((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))
			ENDFOREACH(MATVER)
		ENDIF((NOT DEFINED MATLAB_ROOT) OR ("${MATLAB_ROOT}" STREQUAL ""))

		FIND_PROGRAM(MATLAB_MEX
			mex
			PATHS
			${MATLAB_ROOT}/bin
		)

	ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		# On a Linux system.  The goal is to find MATLAB_ROOT.
		SET(MATLAB_MEXEXT .mexa)
		SET(LIBRARY_EXTENSION .so)

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
		
			FIND_PROGRAM(MATLAB_MEX_POSSIBLE_LINK
			NAMES
			mex
			PATHS
			${MATLAB_ROOT}/bin
			/opt/MATLAB/R2019a/bin
			/opt/MATLAB/R2018b/bin
			/opt/MATLAB/R2018a/bin
			/opt/MATLAB/R2017b/bin
			/opt/MATLAB/R2017a/bin
			/opt/MATLAB/R2016b/bin
			/opt/MATLAB/R2016a/bin
			/opt/MATLAB/R2015b/bin
			/opt/MATLAB/R2015a/bin
			/opt/MATLAB/R2014b/bin
			/opt/MATLAB/R2014a/bin
			/opt/MATLAB/R2013b/bin
			/opt/MATLAB/R2013a/bin
			/opt/MATLAB/R2012b/bin
			/opt/MATLAB/R2012a/bin
			/opt/matlab/bin
			/usr/local/MATLAB/R2019a/bin
			/usr/local/MATLAB/R2018b/bin
			/usr/local/MATLAB/R2018a/bin
			/usr/local/MATLAB/R2017b/bin
			/usr/local/MATLAB/R2017a/bin
			/usr/local/MATLAB/R2016b/bin
			/usr/local/MATLAB/R2016a/bin
			/usr/local/MATLAB/R2015b/bin
			/usr/local/MATLAB/R2015a/bin
			/usr/local/MATLAB/R2014b/bin
			/usr/local/MATLAB/R2014a/bin
			/usr/local/MATLAB/R2013b/bin
			/usr/local/MATLAB/R2013a/bin
			/usr/local/MATLAB/R2012b/bin
			/usr/local/MATLAB/R2012a/bin
			/usr/local/matlab/bin
			$ENV{HOME}/MATLAB/R2019a/bin
			$ENV{HOME}/MATLAB/R2018b/bin
			$ENV{HOME}/MATLAB/R2018a/bin
			$ENV{HOME}/MATLAB/R2017b/bin
			$ENV{HOME}/MATLAB/R2017a/bin
			$ENV{HOME}/MATLAB/R2016b/bin
			$ENV{HOME}/MATLAB/R2016a/bin
			$ENV{HOME}/MATLAB/R2015b/bin
			$ENV{HOME}/MATLAB/R2015a/bin
			$ENV{HOME}/MATLAB/R2014b/bin
			$ENV{HOME}/MATLAB/R2014a/bin
			$ENV{HOME}/MATLAB/R2013b/bin
			$ENV{HOME}/MATLAB/R2013a/bin
			$ENV{HOME}/MATLAB/R2012b/bin
			$ENV{HOME}/MATLAB/R2012a/bin
			$ENV{HOME}/matlab/bin
			NO_DEFAULT_PATH
			)
			# Check if the mex link is found
			if(NOT "${MATLAB_MEX_POSSIBLE_LINK}" STREQUAL "MATLAB_MEX_POSSIBLE_LINK-NOTFOUND" AND NOT "${MATLAB_MEX_POSSIBLE_LINK}" STREQUAL "")
				get_filename_component(MATLAB_MEX "${MATLAB_MEX_POSSIBLE_LINK}" REALPATH)
				get_filename_component(MATLAB_BIN_ROOT "${MATLAB_MEX}" PATH)
				# Strip ./bin/.
				GET_FILENAME_COMPONENT(MATLAB_ROOT "${MATLAB_BIN_ROOT}" PATH)
			else()
				set(MATLAB_ROOT "")
				set(MATLAB_MEX "")
			endif()
		endif(NOT ${locate_matlab_result} AND NOT "${locate_matlab_output}" STREQUAL "")
		else()
			if(EXISTS  "${MATLAB_ROOT}/bin/mex")
				set(MATLAB_MEX "${MATLAB_ROOT}/bin/mex")
			else(EXISTS  "${MATLAB_ROOT}/bin/mex")
				set(MATLAB_ROOT "")
				set(MATLAB_MEX "")
			endif(EXISTS  "${MATLAB_ROOT}/bin/mex")
		endif()
	ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
ENDIF(WIN32)

IF("${MATLAB_MEX}" STREQUAL "")
	MESSAGE(STATUS "Could not find MATLAB mex compiler; try specifying NIHU_MATLAB_PATH if Matlab is already installed.")
ELSE("${MATLAB_MEX}" STREQUAL "")
	# Check if system is 64 bit
	if(CMAKE_SIZEOF_VOID_P EQUAL 8) 
		set(MATLAB_MEXEXT "${MATLAB_MEXEXT}64")
	# Otherwise 32 bit
	else(CMAKE_SIZEOF_VOID_P EQUAL 8)
		set(MATLAB_MEXEXT "${MATLAB_MEXEXT}32")
	endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

	MESSAGE(STATUS "Found MATLAB mex compiler: ${MATLAB_MEX}")
	MESSAGE(STATUS "MATLAB root: ${MATLAB_ROOT}")
	MESSAGE(STATUS "MEX file extension is: ${MATLAB_MEXEXT}")
	SET(MATLABMEX_FOUND 1)
ENDIF("${MATLAB_MEX}" STREQUAL "")

# MARK_AS_ADVANCED(
# 	MATLAB_MEX
# 	MATLABMEX_FOUND
# 	MATLAB_ROOT
# )
