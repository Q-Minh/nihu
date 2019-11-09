# NiHuTesting.cmake

message(STATUS "Configuring testing ...")

if (NIHU_BUILD_GTEST_TESTS)
	message(STATUS "\tGTest test targets enabled")
	add_subdirectory(gtest)
endif()

# Check if install 
if(NIHU_ENABLE_TEST_INSTALL AND NOT DEFINED NIHU_DISABLE_INSTALL)
	if(${NIHU_INSTALL_DIR} MATCHES ${CMAKE_BINARY_DIR})
		set(NIHU_ENABLE_TEST_INSTALL 0)
		message(STATUS "\tInstallation of tests disabled (build dir is same as install dir)")
	else()
		message(STATUS "\tInstallation of tests enabled")
	endif()
endif()
# Enable testing globally
enable_testing()
# Check if automatic build of tests are disabled
if(NOT NIHU_DISABLE_TEST_BUILD)
	add_subdirectory(test)
else()
	add_subdirectory(test EXCLUDE_FROM_ALL)
endif()

message(STATUS "\tTesting configured successfully")
