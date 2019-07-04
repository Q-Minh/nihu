# cmake FFTW for NiHu

# Options:
#	NIHU_FFTW_INSTALL
#		When set to non-zero 
# 	NIHU_FFTW_PATH 		
#		Specify a path where to look for FFTW headers and libraries
#	NIHU_FFTW_ARCHIVE
#		Specify an archive for extracting FFTW
#	NIHU_FFTW_VERSION	
#		Specify which FFTW version to use, e.g. "3.3.5"
#		Currently supported versions: 3.3.5

# Download location: ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll64.zip

# Supported FFTW versions
set(NIHU_SUPPORTED_FFTW_VERSIONS "3.3.5")

# The default Eigen version
set(NIHU_DEFAULT_FFTW_VERSION "3.3.5")

# Process fftw path
if (NOT DEFINED NIHU_FFTW_INSTALL)
	if (DEFINED NIHU_FFTW_PATH)
		set(FFTW3_DIRECTORY ${NIHU_FFTW_PATH})
	endif ()

	find_package(FFTW)
	if (FFTW3_FOUND)
		message(STATUS "fftw3 was found, include directory: ${FFTW3_INCLUDE_DIRS}")
	endif ()
endif ()

# If FFTW not found, install
if (NOT FFTW3_FOUND OR DEFINED NIHU_FFTW_INSTALL)
	if (NOT WIN32)
		message(FATAL_ERROR "FFTW automatic install only supported for windows")
	endif ()
	
	# Supported FFTW versions
	if(DEFINED NIHU_FFTW_VERSION)
		# Check if the defined FFTW version is supported
		list (FIND NIHU_SUPPORTED_FFTW_VERSIONS ${NIHU_FFTW_VERSION} _index)
		if (${_index} GREATER -1)
			message(STATUS "Selected FFTW version ${NIHU_FFTW_VERSION} is supported")
		else(${_index} GREATER -1)
			message(FATAL_ERROR "FFTW version ${NIHU_FFTW_VERSION} is unsupported, supported versions: ${NIHU_SUPPORTED_FFTW_VERSIONS}")
		endif(${_index} GREATER -1)
	else(DEFINED NIHU_FFTW_VERSION)
		# If undefined, set to default
		set(NIHU_FFTW_VERSION "${NIHU_DEFAULT_FFTW_VERSION}")
	endif(DEFINED NIHU_FFTW_VERSION)

	# Configure FFTW directories
	set(FFTW_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/fftw-${NIHU_FFTW_VERSION}")
	set(FFTW_INSTALL_DIR "${CMAKE_BINARY_DIR}/lib/fftw-${NIHU_FFTW_VERSION}")
	
	# md5 checksum of the downloaded 
	if (NIHU_FFTW_VERSION STREQUAL "3.3.5")
		if (NIHU_SYS_BITS STREQUAL "32")
			set(FFTW_MD5 "f9928318c8a35fa8b594aa75af9f0550") #3.3.5 32bit
		endif ()
		if (NIHU_SYS_BITS STREQUAL "64")
			set(FFTW_MD5 "cb3c5ad19a89864f036e7a2dd5be168c") #3.3.5 64bit
		endif ()
	endif (NIHU_FFTW_VERSION STREQUAL "3.3.5")
	
	# Check if NIHU_EIGEN_TARBALL is set
	if (NOT DEFINED NIHU_FFTW_ARCHIVE)
		# FFTW download path
		set(FFTW_URL "ftp://ftp.fftw.org/pub/fftw/fftw-${NIHU_FFTW_VERSION}-dll${NIHU_SYS_BITS}.zip")
		message(STATUS "FFTW ${NIHU_FFTW_VERSION} headers and libraries will be installed as a part of NiHu")
		set(FFTW_DL_FILE "${CMAKE_SOURCE_DIR}/ThirdParty/fftw-${NIHU_FFTW_VERSION}.zip")
		
		# Check if downloaded file exists and has the correct MD5SUM
		if (EXISTS "${FFTW_DL_FILE}")
			file(MD5 "${FFTW_DL_FILE}" FFTW_EXISTING_MD5)
			if (FFTW_EXISTING_MD5 STREQUAL "${FFTW_MD5}")
				set(FFTW_DOWNLOAD 0)
			else ()
				message(WARNING "MD5 sum of existing fftw3 archive does not match the expected")
				set(FFTW_DOWNLOAD 1)
			endif()
		else()
			set(FFTW_DOWNLOAD 1)
		endif()
		
		# Download as needed
		if (FFTW_DOWNLOAD)
			message(STATUS "Downloading FFTW3 from ${FFTW_URL}")
			file(
				DOWNLOAD "${FFTW_URL}" "${FFTW_DL_FILE}"
				EXPECTED_MD5 "${FFTW_MD5}"
				STATUS FFTW_DL_STATUS
			)
			
			# Check download status
			LIST(GET ${FFTW_DL_STATUS} 1 FFTW_DL_ERROR)
			if(FFTW_DL_ERROR)
				message(FATAL_ERROR "FFTW download failed, cmake will exit")
			endif(FFTW_DL_ERROR)
		endif ()
			message(STATUS "Using existing FFTW archive ${FFTW_DL_FILE}")
	else ()
		# Install eigen from a predownloaded tarball file
		if (EXISTS ${NIHU_FFTW_ARCHIVE})
			set(FFTW_DL_FILE ${NIHU_FFTW_ARCHIVE})
		else ()
			message(FATAL_ERROR "Given FFTW tarball ${NIHU_FFTW_ARCHIVE} does not exist, cmake will exit")
		endif ()
	endif ()
	
	# We have a file to extract
	message(STATUS "Extracting FFTW3 into ${FFTW_SOURCE_DIR}")
	# Create directory for FFTW
	file(MAKE_DIRECTORY ${FFTW_SOURCE_DIR})
	# Unpack FFTW
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xvf "${FFTW_DL_FILE}"
		WORKING_DIRECTORY "${FFTW_SOURCE_DIR}"
		OUTPUT_QUIET)
		
	# Create targets for generating the .lib files 
	if (MINGW)
		add_custom_target("libfftw3-3"
			COMMAND dlltool -d libfftw3-3.def -l libfftw3-3.lib  
			COMMAND dlltool -d libfftw3f-3.def -l libfftw3f-3.lib  
			COMMAND dlltool -d libfftw3l-3.def -l libfftw3l-3.lib  
			WORKING_DIRECTORY "${FFTW_SOURCE_DIR}"
		)
	elseif (MSVC)
		add_custom_target("libfftw3-3" 
			COMMAND lib /def:libfftw3-3.def /machine:X${NIHU_SYS_BITS} /out:fftw3-3.lib
			COMMAND lib /def:libfftw3f-3.def /machine:X${NIHU_SYS_BITS} /out:fftw3f-3.lib
			COMMAND lib /def:libfftw3l-3.def /machine:X${NIHU_SYS_BITS} /out:fftw3l-3.lib
			WORKING_DIRECTORY "${FFTW_SOURCE_DIR}"
		)
	endif ()
	
	# Find the package
	set(FFTW3_DIRECTORY ${FFTW_SOURCE_DIR})
	find_package(FFTW)
	
	if (NOT FFTW3_FOUND)
		message(FATAL_ERROR "Could not find the new FFTW installation")
	endif ()
endif ()

# At this point FFTW should be successfully updated
# Add FFTW include directory
include_directories(${FFTW3_INCLUDE_DIRS})
