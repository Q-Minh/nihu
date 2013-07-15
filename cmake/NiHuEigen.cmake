# Check for eigen
find_package (Eigen REQUIRED)

# Check if eigen found
if(EIGEN_FOUND)
	include_directories ("${EIGEN_INCLUDE_DIRS}")
else(EIGEN_FOUND)
	# Try to install eigen
	message(STATUS "Eigen3 not installed yet, adding as a separate project" )
	# Eigen 3.1.2 download location
	# The following are specific for eigen 3.1.2
	set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/3.1.2.tar.bz2")
	# md5 checksum of the downloaded file tar.bz2
	set(EIGEN_MD5 "e9c081360dde5e7dcb8eba3c8430fde2") 

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
	include_directories("${EIGEN_INCLUDE_DIRS}")

endif(EIGEN_FOUND)