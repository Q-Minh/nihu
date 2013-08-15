# If Eigen headers are set
# (This is the default behaviour on windows.)

# Check for eigen
find_package (Eigen)

# Check if eigen found
if(NOT EIGEN_FOUND)
	
	# Eigen download path
	set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/3.1.2.tar.bz2")
	# md5 checksum of the downloaded file tar.bz2
	set(EIGEN_MD5 "e9c081360dde5e7dcb8eba3c8430fde2") #3.1.2

	# Install as headers	
	if(NOT NIHU_EIGEN_AS_PROJECT)
		message(STATUS "Eigen3 headers will be installed as a part of NiHu")

		set(EIGEN_DL_FILE "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-3.1.2.tar.bz2")
		set(EIGEN_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/eigen-3.1.2")
		set(EIGEN_TEMP_DIR "eigen-eigen-5097c01bcdc4")
		set(EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/include/eigen-3.1.2")

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
			
		# Install the header files
		message(STATUS "Installing Eigen3 headers into ${EIGEN_INSTALL_DIR}")
		execute_process(
			COMMAND ${CMAKE_COMMAND} -E make_directory "${EIGEN_INSTALL_DIR}"
			OUTPUT_QUIET)
		file(COPY "${EIGEN_SOURCE_DIR}/Eigen" DESTINATION "${EIGEN_INSTALL_DIR}/")
		
		# Install the header files
		message(STATUS "Eigen3 headers successfully installed")
		
		# Add the include directory
		set(EIGEN_INCLUDE_DIRS "${EIGEN_INSTALL_DIR}")

	# Otherwise, copy 
	else(NOT NIHU_EIGEN_AS_PROJECT)

		# Try to install eigen
		message(STATUS "Eigen3 not installed yet, will be compiled from source" )

		# Setup source and build directories
		set(EIGEN_SOURCE_PATH "${NIHU_THIRDPARTY_DIR}/eigen-3.1.2")
		set(EIGEN_BUILD_DIR "${CMAKE_BINARY_DIR}/ThirdParty/eigen-3.1.2")
		set(EIGEN_INCLUDES_DIR "${CMAKE_BINARY_DIR}/include/eigen")

		# Include external project library
		include(ExternalProject)
		
		# Add eigen 3 as an external project
		ExternalProject_Add(
			"eigen-3.1.2"
			# Download step
			CMAKE_ARGS 
				"-DCMAKE_INSTALL_PREFIX=${EIGEN_BUILD_DIR}"
				"-DEIGEN_INCLUDE_INSTALL_DIR=${EIGEN_INCLUDES_DIR}"
				"${EIGEN_ADD_CMAKE_ARGS}"
			PREFIX "${CMAKE_BINARY_DIR}/ThirdParty/eigen-3.1.2"
			DOWNLOAD_DIR "${NIHU_THIRDPARTY_DIR}"
			URL "${EIGEN_URL}"
			URL_MD5 "${EIGEN_MD5}"
			# Configuraation step
			SOURCE_DIR "${EIGEN_SOURCE_PATH}"
			BINARY_DIR "${EIGEN_BUILD_DIR}"
			# Logging
			LOG_DOWNLOAD 1
			LOG_CONFIGURE 1
			LOG_BUILD 1
			LOG_INSTALL 1
		)

		# include the install directory of eigen includes
		set(EIGEN_INCLUDE_DIRS "${EIGEN_INCLUDES_DIR}")
	endif(NOT NIHU_EIGEN_AS_PROJECT)

endif(NOT EIGEN_FOUND)
