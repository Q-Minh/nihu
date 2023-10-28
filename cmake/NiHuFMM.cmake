# NiHuFMM.cmake 
# Settings for the FMM module

# NIHU_FMM_DISABLE_PARALLEL
# NIHU_FMM_TRAVERSE 

message(STATUS "Configuring FMM module ...")

# set parallel
if(NOT DEFINED NIHU_FMM_DISABLE_PARALLEL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNIHU_FMM_PARALLEL")
endif() 

# set traversing parameter to fmm_matrix
if(NOT DEFINED NIHU_FMM_TRAVERSE)
	set(NIHU_FMM_TRAVERSE "BFS")
else()
	if((NOT NIHU_FMM_TRAVERSE STREQUAL "BFS")
		AND (NOT NIHU_FMM_TRAVERSE STREQUAL "DFS"))
		message(FATAL_ERROR "\tTraverse parameter NIHU_FMM_TRAVERSE must either be \"BFS\" or \"DFS\"")
	endif()
endif() 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNIHU_FMM_TRAVERSE_${NIHU_FMM_TRAVERSE}")

message(STATUS "\tFMM module configured successfully")
