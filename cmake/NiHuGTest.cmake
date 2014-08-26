# NiHuGTest.cmake
# CMake file for integration of GTest into NiHu

# Check if user specified the directory for GTest
if (DEFINED NIHU_GTEST_PATH)
	# Check if defined directory exists
	if (NOT EXISTS ${NIHU_GTEST_PATH}/include/gtest)
		message(FATAL_ERROR "GTest cannot be found at the specified path: ${NIHU_GTEST_PATH}")
	endif (NOT EXISTS ${NIHU_GTEST_PATH}/include/gtest)
	# GTest has been found successfully
	set (GTEST_FOUND 1)
	# Add the specified path to the include dir list
	set (GTEST_INCLUDE_DIRS NIHU_GTEST_PATH)
else (DEFINED NIHU_GTEST_PATH)
	find_package(GTest)
endif(DEFINED NIHU_GTEST_PATH)

# Check if GTest has been found
if (NOT GTEST_FOUND OR NIHU_GTEST_INSTALL)
	# Configure GTest directories
	set(GTEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/gtest-1.7.0")
	
	# Check if NIHU_EIGEN_TARBALL is set
	if (NOT DEFINED NIHU_GTEST_ARCHIVE)
		# Eigen download path
		set(GTEST_URL "http://developer.nrel.gov/downloads/buildings/openstudio/src/gtest-1.7.0.tar.gz")
		# md5 checksum of the downloaded file
		set(GTEST_MD5 "bad74f626d18f724cad4c9fe2e6ef27d") #1.7.0

		message(STATUS "GTest 1.7.0 headers will be installed as a part of NiHu")

		set(GTEST_DL_FILE "${CMAKE_SOURCE_DIR}/ThirdParty/gtest-1.7.0.tar.gz")
		
		# Download
		message(STATUS "Downloading GTest 1.7.0 from ${GTEST_URL}")
		file(
			DOWNLOAD "${GTEST_URL}" "${GTEST_DL_FILE}"
			EXPECTED_MD5 "${GTEST_MD5}"
			STATUS GTEST_DL_STATUS
		)
		
		# Check download status
		LIST(GET ${GTEST_DL_STATUS} 1 GTEST_DL_ERROR)
		if(EIGEN_DL_ERROR)
			message(FATAL_ERROR "GTest download failed, cmake will exit")
		endif(EIGEN_DL_ERROR)
	else (NOT DEFINED NIHU_GTEST_ARCHIVE)
		# Install gtest from a predownloaded archive file
		if (EXISTS ${NIHU_GTEST_ARCHIVE})
			set(GTEST_DL_FILE ${NIHU_GTEST_ARCHIVE})
		else (EXISTS ${NIHU_GTEST_ARCHIVE})
			message(FATAL_ERROR "Given GTest archive ${NIHU_GTEST_ARCHIVE} does not exist, cmake will exit")
		endif (EXISTS ${NIHU_GTEST_ARCHIVE})
	endif (NOT DEFINED NIHU_GTEST_ARCHIVE)
	
	# Extract
	message(STATUS "Extracting GTest 1.7. into ${GTEST_SOURCE_DIR}")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xvf "${GTEST_DL_FILE}" 
		WORKING_DIRECTORY "${NIHU_THIRDPARTY_DIR}"
		OUTPUT_QUIET)
	
	# Add the include directory
	set(GTEST_MAIN_DIR "${GTEST_SOURCE_DIR}")
	set(GTEST_INCLUDE_DIRS "${GTEST_SOURCE_DIR}/include")
	
endif (NOT GTEST_FOUND OR NIHU_GTEST_INSTALL)