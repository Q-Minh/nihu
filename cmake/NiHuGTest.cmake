# NiHuGTest.cmake
# CMake file for integration of GTest into NiHu

#
# NIHU_GTEST_ARCHIVE
# NIHU_GTEST_INSTALL 
# NIHU_GTEST_PATH

set(NIHU_GTEST_SUPPORTED_VERSIONS "1.7.0")

set(NIHU_GTEST_DEFAULT_VERSION "1.7.0")

message(STATUS "Looking for GTest ...")

if(NOT DEFINED NIHU_GTEST_VERSION)
	set(NIHU_GTEST_VERSION ${NIHU_GTEST_DEFAULT_VERSION})
endif()

# Check if user specified the directory for GTest
if(DEFINED NIHU_GTEST_PATH)
	message(STATUS "\tGTest path (user input): ${NIHU_GTEST_PATH}")
	# Check if defined directory exists
	if(NOT EXISTS ${NIHU_GTEST_PATH}/include/gtest)
		message(FATAL_ERROR "\tGTest cannot be found at the specified path: ${NIHU_GTEST_PATH}")
	endif()
	# GTest has been found successfully
	set(GTEST_FOUND 1)
	# Add the specified path to the include dir list
	set(GTEST_INCLUDE_DIRS NIHU_GTEST_PATH)
else()
	find_package(GTest QUIET)
endif()

# Check if GTest has been found
if(NOT GTEST_FOUND OR NIHU_GTEST_INSTALL)
	# Configure GTest directories
	set(GTEST_SOURCE_DIR "${NIHU_THIRDPARTY_DIR}/gtest-${NIHU_GTEST_VERSION}")
	
	# Check if NIHU_GTEST_ARCHIVE is set
	if(NOT DEFINED NIHU_GTEST_ARCHIVE)
		# Eigen download path
		set(GTEST_URL "http://developer.nrel.gov/downloads/buildings/openstudio/src/gtest-${NIHU_GTEST_VERSION}.tar.gz")
		# md5 checksum of the downloaded file
		if(NIHU_GTEST_VERSION STREQUAL "1.7.0")
			set(GTEST_MD5 "bad74f626d18f724cad4c9fe2e6ef27d") #1.7.0
		endif()

		message(STATUS "\tGTest ${NIHU_GTEST_VERSION} headers will be installed")

		set(GTEST_DL_FILE "${NIHU_THIRDPARTY_DIR}/gtest-${NIHU_GTEST_VERSION}.tar.gz")
		
		# Download
		message(STATUS "\tDownloading GTest ${NIHU_GTEST_VERSION} from ${GTEST_URL}")
		file(
			DOWNLOAD "${GTEST_URL}" "${GTEST_DL_FILE}"
			EXPECTED_MD5 "${GTEST_MD5}"
			STATUS GTEST_DL_STATUS
		)
		
		# Check download status
		list(GET ${GTEST_DL_STATUS} 1 GTEST_DL_ERROR)
		if(GTEST_DL_ERROR)
			message(FATAL_ERROR "\tGTest download failed")
		endif()
	else()
		# Install gtest from a predownloaded archive file
		if(EXISTS ${NIHU_GTEST_ARCHIVE})
			set(GTEST_DL_FILE ${NIHU_GTEST_ARCHIVE})
		else()
			message(FATAL_ERROR "\tThe given GTest archive ${NIHU_GTEST_ARCHIVE} does not exist")
		endif()
	endif()
	
	# Extract
	message(STATUS "\tExtracting GTest ${NIHU_GTEST_VERSION} into ${GTEST_SOURCE_DIR}")
	execute_process(
		COMMAND ${CMAKE_COMMAND} -E tar xvf "${GTEST_DL_FILE}" 
		WORKING_DIRECTORY "${NIHU_THIRDPARTY_DIR}"
		OUTPUT_QUIET)
	
	# Add the include directory
	set(GTEST_MAIN_DIR "${GTEST_SOURCE_DIR}")
	set(GTEST_INCLUDE_DIRS "${GTEST_SOURCE_DIR}/include")
endif()

# Check if gtest directories are found
if(NOT DEFINED GTEST_INCLUDE_DIRS)
	message(STATUS "\tGTest not found, or could not be installed")
	message(STATUS "\tBuild of GTest tests is disabled")
	set(NIHU_BUILD_GTEST_TESTS 0)
else ()
	message(STATUS "\tGTest headers will be imported from: ${GTEST_INCLUDE_DIRS}")
	message(STATUS "\tBuild of GTest test is enabled")
	set(NIHU_BUILD_GTEST_TESTS 1)
endif()
