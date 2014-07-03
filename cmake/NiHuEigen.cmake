# If Eigen headers are set
# (This is the default behaviour on windows.)

# Check for given eigen path
if(DEFINED NIHU_EIGEN_PATH)
	set(EIGEN_FOUND 1)
	# Check if Eigen directory exists
	if(NOT EXISTS "${NIHU_EIGEN_PATH}/Eigen")
		message(ERROR "Eigen can not be found at the specified path: ${NIHU_EIGEN_PATH}")
	endif(NOT EXISTS "${NIHU_EIGEN_PATH}/Eigen")
	set(EIGEN_INCLUDE_DIRS NIHU_EIGEN_PATH)
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
	
	# Eigen download path
	set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/3.2.0.tar.bz2")
	# md5 checksum of the downloaded file tar.bz2
	set(EIGEN_MD5 "894381be5be65bb7099c6fd91d61b357") #3.2.0

	message(STATUS "Eigen3 headers will be installed as a part of NiHu")

	set(EIGEN_DL_FILE "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-3.2.0.tar.bz2")
	set(EIGEN_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-3.2.0")
	set(EIGEN_TEMP_DIR "eigen-eigen-ffa86ffb5570")
	set(EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/include/eigen-3.2.0")

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

#	# Apply eigen patch
#	if(NOT NIHU_EIGEN_DISABLE_PATCH)
#		message(STATUS "Applying Eigen patch")
#		execute_process(COMMAND "patch" "--directory=${EIGEN_SOURCE_DIR}" "--strip=2" "--input=${NIHU_THIRDPARTY_DIR}/patches/eigen-3.2.0.patch")
#	endif(NOT NIHU_EIGEN_DISABLE_PATCH)
			
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
