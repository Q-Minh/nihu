cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(NiHu)

set(NIHU_VERSION_MAJOR 2)
set(NIHU_VERSION_MINOR 0)

option(NIHU_ENABLE_MATLAB "Enable MATLAB" OFF)
option(NIHU_ENABLE_TESTING "Include unit tests" OFF)
option(NIHU_ENABLE_FMM "Enable fast multipole method" ON)
option(NIHU_FMM_DISABLE_PARALLEL "Disable parallel Fast Multipole Method" OFF)
option(NIHU_ENABLE_TUTORIALS "Enable NiHu tutorials" OFF)
option(NIHU_ENABLE_APPLICATIONS "Enable NiHu applications" OFF)
option(NIHU_ENABLE_DOC "Enable NiHu documentation" OFF)

include(cmake/NiHuSystem.cmake)

if (NOT TARGET Boost::boost)
    find_package(Boost REQUIRED)
endif()

if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3...3.4 CONFIG REQUIRED)
endif()

if (NOT TARGET FFTW3::fftw3)
    find_package(FFTW3 CONFIG REQUIRED)
endif()

if (NOT TARGET OpenMP::OpenMP_CXX)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

add_subdirectory(nihu)

include(${CMAKE_SOURCE_DIR}/cmake/NiHuInstall.cmake)
