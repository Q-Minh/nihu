# message(STATUS "Configuring installation steps ...")

# ### Setup directories for installation
# if(NOT DEFINED NIHU_INSTALL_DIR)
# 	set(NIHU_INSTALL_DIR ${CMAKE_BINARY_DIR})
# endif(NOT DEFINED NIHU_INSTALL_DIR)

# set(CMAKE_INSTALL_PREFIX ${NIHU_INSTALL_DIR})

# # Select all hpp files for installation
# set(NIHU_HPP_DIRECTORIES 
# 	"core" ; "interface" ; "library" ; "tmp" ; "util" ; "aca")
# if (NOT DEFINED NIHU_DISABLE_FMM)
# 	set(NIHU_HPP_DIRECTORIES "${NIHU_HPP_DIRECTORIES}" ; "fmm")
# endif()

# foreach(HPP_DIRECTORY ${NIHU_HPP_DIRECTORIES})
# 	install(DIRECTORY ${HPP_DIRECTORY} DESTINATION include FILES_MATCHING PATTERN "*.hpp")
# endforeach(HPP_DIRECTORY)

# # Matlab installation section
# set(NIHU_MATLAB_DIRECTORIES
# 	"analytic" ; "compatibility" ; "meshing")

# foreach(MATLAB_DIRECTORY ${NIHU_MATLAB_DIRECTORIES})
# 	install(DIRECTORY "matlab/${MATLAB_DIRECTORY}" DESTINATION matlab FILES_MATCHING PATTERN "*.m")
# endforeach(MATLAB_DIRECTORY)

# # Installation rule for matlab demos
# install(DIRECTORY "matlab/nihudemos" DESTINATION matlab)

# # Create installation rule for matlab install script
# install(FILES "matlab/install.m" DESTINATION matlab)

# # documentation installation
# if(NOT DEFINED NIHU_DISABLE_DOC)
# 	if(NOT (${NIHU_INSTALL_PATH} MATCHES ${CMAKE_BINARY_DIR}))
# 		message(STATUS "\tDocumentation will be installed into ${CMAKE_INSTALL_PREFIX}/doc")
# 		install(DIRECTORY ${NIHU_HTML_DOC_DIR} DESTINATION doc)
# 	endif()
# endif()

# message(STATUS "\tInstallation steps configured successfully")
