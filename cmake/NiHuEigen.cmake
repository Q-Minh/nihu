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
set(NIHU_EIGEN_SUPPORTED_VERSIONS "3.2.7" "3.2.9" "3.3.7")

# The default Eigen version
set(NIHU_EIGEN_DEFAULT_VERSION "3.3.7")

message(STATUS "Looking for Eigen3 ...")

# Check for given eigen path
if(DEFINED NIHU_EIGEN_PATH)
    message(STATUS "\tEigen path (user input): ${NIHU_EIGEN_PATH}")
    # Check if Eigen directory exists
    if(NOT EXISTS "${NIHU_EIGEN_PATH}/Eigen")
        message(FATAL_ERROR "Eigen can not be found at the specified path: ${NIHU_EIGEN_PATH}")
    endif()
    # Eigen found successfully 
    set(EIGEN_FOUND 1)
    set(EIGEN_INCLUDE_DIRS "${NIHU_EIGEN_PATH}")
# If the path is not defined
else()
    if(UNIX)
        if(NOT NIHU_EIGEN_INSTALL)
            find_package(Eigen)
        endif()
    else()
        if(WIN32)
            set(EIGEN_FOUND 0)
            set(NIHU_EIGEN_INSTALL 1)
        endif()
    endif()
endif()

# Check if eigen found
if(NOT EIGEN_FOUND OR NIHU_EIGEN_INSTALL)
    # Supported Eigen versions
    if(DEFINED NIHU_EIGEN_VERSION)
        # Check if the defined Eigen version is supported
        list (FIND NIHU_EIGEN_SUPPORTED_VERSIONS ${NIHU_EIGEN_VERSION} _index)
        if (${_index} GREATER -1)
            message(STATUS "\tSelected Eigen version ${NIHU_EIGEN_VERSION} is supported")
        else()
            message(FATAL_ERROR "\tEigen version ${NIHU_EIGEN_VERSION} is unsupported, supported versions: ${NIHU_EIGEN_SUPPORTED_VERSIONS}")
        endif()
    else()
        set(NIHU_EIGEN_VERSION "${NIHU_EIGEN_DEFAULT_VERSION}")
    endif()

    # Configure Eigen directories
    set(EIGEN_SOURCE_DIR "${NIHU_THIRDPARTY_DIR}/eigen-${NIHU_EIGEN_VERSION}")
    if(NIHU_EIGEN_VERSION STREQUAL "3.2.7")
        set(EIGEN_TEMP_DIR "eigen-eigen-b30b87236a1b")
        set(EIGEN_MD5 "cc1bacbad97558b97da6b77c9644f184") #3.2.7
    elseif(NIHU_EIGEN_VERSION STREQUAL "3.2.9")
        set(EIGEN_TEMP_DIR "eigen-eigen-dc6cfdf9bcec")
        set(EIGEN_MD5 "de11bfbfe2fd2dc4b32e8f416f58ee98") #3.2.9
    elseif(NIHU_EIGEN_VERSION STREQUAL "3.3.7")
        set(EIGEN_TEMP_DIR "eigen-eigen-323c052e1731")
        set(EIGEN_MD5 "05b1f7511c93980c385ebe11bd3c93fa") #3.3.7
    endif()
    
    set(EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/include/eigen-${NIHU_EIGEN_VERSION}")
    
    # Check if NIHU_EIGEN_ARCHIVE is set
    if(NOT DEFINED NIHU_EIGEN_ARCHIVE)
        # Eigen download path
        set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/${NIHU_EIGEN_VERSION}.tar.bz2")
        # md5 checksum of the downloaded file tar.bz2
        message(STATUS "\tEigen ${NIHU_EIGEN_VERSION} headers will be installed")

        set(EIGEN_DL_FILE "${NIHU_THIRDPARTY_DIR}/eigen-${NIHU_EIGEN_VERSION}.tar.bz2")
        
        # Check if downloaded file exists
        if(EXISTS "${EIGEN_DL_FILE}")
            file(MD5 "${EIGEN_DL_FILE}" EIGEN_EXISTING_MD5)
            if (EIGEN_EXISTING_MD5 STREQUAL "${EIGEN_MD5}")
                set(EIGEN_DOWNLOAD 0)
            else ()
                message(WARNING "\tMD5 sum of existing Eigen archive does not match expected")
                set(EIGEN_DOWNLOAD 1)
            endif()
        else()
            set(EIGEN_DOWNLOAD 1)
        endif()
        
        # Download as needed
        if(EIGEN_DOWNLOAD)
            message(STATUS "\tDownloading Eigen3 from ${EIGEN_URL}")
            file(
                DOWNLOAD "${EIGEN_URL}" "${EIGEN_DL_FILE}"
                EXPECTED_MD5 "${EIGEN_MD5}"
                STATUS EIGEN_DL_STATUS
            )
            
            # Check download status
            list(GET ${EIGEN_DL_STATUS} 1 EIGEN_DL_ERROR)
            if(EIGEN_DL_ERROR)
                message(FATAL_ERROR "\tEigen download failed")
            endif(EIGEN_DL_ERROR)
        else()
            message(STATUS "\tUsing existing Eigen archive ${EIGEN_DL_FILE}")
        endif()
    else()
        # Install eigen from a predownloaded tarball file
        if(EXISTS ${NIHU_EIGEN_TARBALL})
            set(EIGEN_DL_FILE ${NIHU_EIGEN_TARBALL})
        else()
            message(FATAL_ERROR "\tGiven Eigen tarball ${NIHU_EIGEN_TARBALL} does not exist")
        endif()
    endif()
    
    # Extract
    message(STATUS "\tExtracting Eigen3 into ${EIGEN_SOURCE_DIR}")
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

    # Copy the header files
    message(STATUS "\tCopying Eigen3 headers into ${EIGEN_INSTALL_DIR}")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E make_directory "${EIGEN_INSTALL_DIR}"
        OUTPUT_QUIET)
    file(COPY "${EIGEN_SOURCE_DIR}/Eigen" DESTINATION "${EIGEN_INSTALL_DIR}/")
        
    # Install the header files
    message(STATUS "\tEigen3 headers successfully copied")
    
    # Add the install rule for Eigen
    if(NOT (${NIHU_INSTALL_DIR} MATCHES ${CMAKE_BINARY_DIR}))
        install(DIRECTORY "${EIGEN_SOURCE_DIR}/Eigen" DESTINATION include)
    endif()

    # Add the include directory
    set(EIGEN_INCLUDE_DIRS "${EIGEN_INSTALL_DIR}")
endif()
