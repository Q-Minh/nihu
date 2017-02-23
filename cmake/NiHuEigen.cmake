# NiHuEigen.cmake
# 
# This file is used for setting up Eigen directories during NiHu install.
# Two alternative ways of using Eigen are possible:
# (1) If you have Eigen headers installed on your computer, then cmake will 
#     autmatically look for the Eigen path, this is the default behavior on Unix
#     systems. If you do not want to use your installed version, you can force
#     cmake to install Eigen as a part of NiHu by setting the option
#     -DNIHU_EIGEN_INSTALL=1
# (2) If you do not have Eigen headers installed, Eigen is automatically 
#     downloaded and installed as a part of NiHu. This is the default behavior
#     on Windows systems.

# Supported Eigen versions
set(NIHU_SUPPORTED_EIGEN_VERSIONS "3.2.0" "3.2.7" "3.2.9")

# The default Eigen version
set(NIHU_DEFAULT_EIGEN_VERSION "3.2.9")

# Check for given eigen path
if(DEFINED NIHU_EIGEN_PATH)
	# Check if Eigen directory exists
	if(NOT EXISTS "${NIHU_EIGEN_PATH}/Eigen")
		message(FATAL_ERROR "Eigen can not be found at the specified path: ${NIHU_EIGEN_PATH}")
	endif(NOT EXISTS "${NIHU_EIGEN_PATH}/Eigen")
	# Eigen found successfully 
	set(EIGEN_FOUND 1)
	set(EIGEN_INCLUDE_DIRS "${NIHU_EIGEN_PATH}")
# If the path is not defined
else(DEFINED NIHU_EIGEN_PATH)
	if(UNIX)
		if(NOT NIHU_EIGEN_INSTALL)
			find_package(Eigen)
		endif(NOT NIHU_EIGEN_INSTALL)
	else(UNIX)
		if(WIN32)
			set(EIGEN_FOUND 0)
			set(NIHU_EIGEN_INSTALL 1)
		endif(WIN32)
	endif(UNIX)
endif(DEFINED NIHU_EIGEN_PATH)

# Check if eigen found
if(NOT EIGEN_FOUND OR NIHU_EIGEN_INSTALL)
	# Supported Eigen versions
	if(DEFINED NIHU_EIGEN_VERSION)
		# Check if the defined Eigen version is supported
		list (FIND NIHU_SUPPORTED_EIGEN_VERSIONS ${NIHU_EIGEN_VERSION} _index)
		if (${_index} GREATER -1)
			message(STATUS "Selected Eigen version ${NIHU_EIGEN_VERSION} is supported")
		else(${_index} GREATER -1)
			message(FATAL_ERROR "Eigen version ${NIHU_EIGEN_VERSION} is unsupported, supported versions: ${NIHU_SUPPORTED_EIGEN_VERSIONS}")
		endif(${_index} GREATER -1)
	else(DEFINED NIHU_EIGEN_VERSION)
		set(NIHU_EIGEN_VERSION "${NIHU_DEFAULT_EIGEN_VERSION}")
	endif(DEFINED NIHU_EIGEN_VERSION)

	# Configure Eigen directories
	set(EIGEN_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-${NIHU_EIGEN_VERSION}")
	if (NIHU_EIGEN_VERSION STREQUAL "3.2.0")
		set(EIGEN_TEMP_DIR "eigen-eigen-ffa86ffb5570")
		set(EIGEN_MD5 "894381be5be65bb7099c6fd91d61b357") #3.2.0
	elseif(NIHU_EIGEN_VERSION STREQUAL "3.2.7")
		set(EIGEN_TEMP_DIR "eigen-eigen-b30b87236a1b")
		set(EIGEN_MD5 "cc1bacbad97558b97da6b77c9644f184") #3.2.7
	elseif(NIHU_EIGEN_VERSION STREQUAL "3.2.9")
		set(EIGEN_TEMP_DIR "eigen-eigen-dc6cfdf9bcec")
		set(EIGEN_MD5 "de11bfbfe2fd2dc4b32e8f416f58ee98") #3.2.9
	endif()
	
	set(EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/include/eigen-${NIHU_EIGEN_VERSION}")
	
	# Check if NIHU_EIGEN_TARBALL is set
	if (NOT DEFINED NIHU_EIGEN_TARBALL)
		# Eigen download path
		set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/${NIHU_EIGEN_VERSION}.tar.bz2")
		# md5 checksum of the downloaded file tar.bz2
		message(STATUS "Eigen ${NIHU_EIGEN_VERSION} headers will be installed as a part of NiHu")

		set(EIGEN_DL_FILE "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-${NIHU_EIGEN_VERSION}.tar.bz2")
		
		# Download
		message(STATUS "Downloading Eigen3 from ${EIGEN_URL}")
		file(
			DOWNLOAD "${EIGEN_URL}" "${EIGEN_DL_FILE}"
			EXPECTED_MD5 "${EIGEN_MD5}"
			STATUS EIGEN_DL_STATUS
		)
		
		# Check download status
		LIST(GET ${EIGEN_DL_STATUS} 1 EIGEN_DL_ERROR)
		if(EIGEN_DL_ERROR)
			message(FATAL_ERROR "Eigen download failed, cmake will exit")
		endif(EIGEN_DL_ERROR)
	else ()
		# Install eigen from a predownloaded tarball file
		if (EXISTS ${NIHU_EIGEN_TARBALL})
			set(EIGEN_DL_FILE ${NIHU_EIGEN_TARBALL})
		else ()
			message(FATAL_ERROR "Given Eigen tarball ${NIHU_EIGEN_TARBALL} does not exist, cmake will exit")
		endif ()
	endif ()
	
	# Extract
	message(STATUS "Extracting Eigen3 into ${EIGEN_SOURCE_DIR}")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xvf "${EIGEN_DL_FILE}" 
		WORKING_DIRECTORY "${NIHU_THIRDPARTY_DIR}"
		OUTPUT_QUIET)
	
	# Rename directory
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E remove_directory "${EIGEN_SOURCE_DIR}"
		OUTPUT_QUIET)
			
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E rename "${NIHU_THIRDPARTY_DIR}/${EIGEN_TEMP_DIR}" "${EIGEN_SOURCE_DIR}"
		OUTPUT_QUIET)

	# Apply eigen patch
#	if(UNIX AND NOT NIHU_EIGEN_DISABLE_PATCH)
#		message(STATUS "Applying Eigen patch")
#		execute_process(COMMAND "patch" "--directory=${EIGEN_SOURCE_DIR}" "--strip=2" "--input=${NIHU_THIRDPARTY_DIR}/patches/eigen-3.2.0.patch")
#	endif(UNIX AND NOT NIHU_EIGEN_DISABLE_PATCH)
			
	# Copy the header files
	message(STATUS "Copying Eigen3 headers into ${EIGEN_INSTALL_DIR}")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E make_directory "${EIGEN_INSTALL_DIR}"
		OUTPUT_QUIET)
	file(COPY "${EIGEN_SOURCE_DIR}/Eigen" DESTINATION "${EIGEN_INSTALL_DIR}/")
		
	# Install the header files
	message(STATUS "Eigen3 headers successfully copied")
	
	# Add the install rule for Eigen
	if(NOT (${NIHU_INSTALL_DIR} MATCHES ${CMAKE_BINARY_DIR}))
		install(DIRECTORY "${EIGEN_SOURCE_DIR}/Eigen" DESTINATION include)
	endif()

	# Add the include directory
	set(EIGEN_INCLUDE_DIRS "${EIGEN_INSTALL_DIR}")
	
endif(NOT EIGEN_FOUND OR NIHU_EIGEN_INSTALL)
