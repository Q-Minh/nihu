include(GNUInstallDirs)

install(
	TARGETS nihu
	EXPORT nihuTargets 
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

set(_nihu_install_folders)
list(APPEND _nihu_install_folders 
	"aca" 
	"core" 
	"cqm" 
	"fmm" 
	"library" 
	"tmp" 
	"util")

foreach(_nihu_install_folder ${_nihu_install_folders})
	install(
		DIRECTORY ${CMAKE_SOURCE_DIR}/nihu/${_nihu_install_folder}
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nihu
		FILES_MATCHING 
			PATTERN "*.h"
			PATTERN "*.hpp"
	)	
endforeach()

# install(
# 	DIRECTORY ${CMAKE_SOURCE_DIR}/nihu
# 	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
# 	PATTERN "applications" EXCLUDE
# 	PATTERN "demo" EXCLUDE
# 	PATTERN "gtest" EXCLUDE
# 	PATTERN "interface" EXCLUDE
# 	PATTERN "matlab" EXCLUDE
# 	PATTERN "matlab_util" EXCLUDE
# 	PATTERN "test" EXCLUDE
# 	PATTERN "tutorial" EXCLUDE
# 	FILES_MATCHING 
# 		PATTERN "*.h"
# 		PATTERN "*.hpp"
# )

install(
	EXPORT nihuTargets
	FILE nihuTargets.cmake
	NAMESPACE nihu::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nihu
	)
	
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/nihuConfigVersion.cmake
    VERSION "${NIHU_VERSION_MAJOR}.${NIHU_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
    )

configure_package_config_file(
	${CMAKE_SOURCE_DIR}/cmake/nihuConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nihuConfig.cmake
INSTALL_DESTINATION 
	${CMAKE_INSTALL_LIBDIR}/cmake/nihu
)

install(
FILES 
	"${CMAKE_CURRENT_BINARY_DIR}/nihuConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/nihuConfigVersion.cmake"
DESTINATION 
	${CMAKE_INSTALL_LIBDIR}/cmake/nihu
	)
