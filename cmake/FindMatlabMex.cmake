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
	# This is untested but taken from the older FindMatlab.cmake script as well as
	# the modifications by Ramon Casero and Tom Doel for Gerardus.

	# Search for a version of Matlab available, starting from the most modern one
	# to older versions.
	FOREACH(MATVER "8.2" "8.1" "8.0" "7.14" "7.13" "7.12" "7.11" "7.10" "7.9" "7.8" "7.7" "7.6" "7.5" "7.4")
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

	message(STATUS "Matlab root is ${MATLAB_ROOT}")
	
	if(EXISTS "${MATLAB_ROOT}/bin/mex.bat")
		set(MATLAB_MEX "${MATLAB_ROOT}/bin/mex.bat")
	endif(EXISTS "${MATLAB_ROOT}/bin/mex.bat")
	
	message(STATUS "Matlab mex is ${MATLAB_MEX}")
ELSE(WIN32)
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
			FOREACH(MATVER "R2013b" "R2013a" "R2012b" "R2012a" "R2011b" "R2011a" "R2010b" "R2010a" "R2009b" "R2009a" "R2008b")
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

		# TODO: REGEXs do not work here somehow
		FIND_PROGRAM(MATLAB_MEX_POSSIBLE_LINK
			NAMES
			mex
			PATHS
			/opt/MATLAB/R2012b/bin
			/opt/MATLAB/R2012a/bin
			${MATLAB_ROOT}/bin
			/opt/matlab/bin
			/usr/local/matlab/bin
			$ENV{HOME}/matlab/bin
			# Now all the versions
			/opt/MATLAB/[rR]20[0-9][0-9][abAB]/bin
			/usr/local/matlab/[rR]20[0-9][0-9][abAB]/bin
			/opt/matlab-[rR]20[0-9][0-9][abAB]/bin
			/opt/matlab_[rR]20[0-9][0-9][abAB]/bin
			/usr/local/matlab-[rR]20[0-9][0-9][abAB]/bin
			/usr/local/matlab_[rR]20[0-9][0-9][abAB]/bin
			$ENV{HOME}/matlab/[rR]20[0-9][0-9][abAB]/bin
			$ENV{HOME}/matlab-[rR]20[0-9][0-9][abAB]/bin
			$ENV{HOME}/matlab_[rR]20[0-9][0-9][abAB]/bin
			NO_DEFAULT_PATH
		)

		# Check if the mex link is found
		IF(NOT ("${MATLAB_MEX_POSSIBLE_LINK}" STREQUAL "MATLAB_MEX_POSSIBLE_LINK-NOTFOUND"))
			GET_FILENAME_COMPONENT(MATLAB_MEX "${MATLAB_MEX_POSSIBLE_LINK}" REALPATH)
			GET_FILENAME_COMPONENT(MATLAB_BIN_ROOT "${MATLAB_MEX}" PATH)
			# Strip ./bin/.
			GET_FILENAME_COMPONENT(MATLAB_ROOT "${MATLAB_BIN_ROOT}" PATH)
		ELSE()
		ENDIF()
	ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
ENDIF(WIN32)

IF("${MATLAB_MEX}" STREQUAL "" AND "${MatlabMex_FIND_REQUIRED}")
	MESSAGE(STATUS "Could not find MATLAB mex compiler; try specifying MATLAB_ROOT.")
ELSE("${MATLAB_MEX}" STREQUAL "" AND "${MatlabMex_FIND_REQUIRED}")
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
ENDIF("${MATLAB_MEX}" STREQUAL "" AND "${MatlabMex_FIND_REQUIRED}")

MARK_AS_ADVANCED(
	MATLAB_MEX
	MATLABMEX_FOUND
	MATLAB_ROOT
)
