# NiHuFMM.cmake 
# Settings for the FMM module

# NIHU_FMM_DISABLE_PARALLEL
# NIHU_FMM_TRAVERSE 

message(STATUS "Configuring FMM module ...")

set(NIHU_FMM_COMPILE_DEFINITIONS "")
# set parallel
if(${NIHU_FMM_DISABLE_PARALLEL})
	list(APPEND NIHU_FMM_COMPILE_DEFINITIONS "-DNIHU_FMM_PARALLEL")
endif() 

# set traversing parameter to fmm_matrix
if((NOT DEFINED NIHU_FMM_TRAVERSE) OR (NIHU_FMM_TRAVERSE STREQUAL "BFS"))
	set(NIHU_FMM_TRAVERSE "BFS")
elseif(
	(NOT NIHU_FMM_TRAVERSE STREQUAL "BFS")
	AND (NOT NIHU_FMM_TRAVERSE STREQUAL "DFS"))
	message(FATAL_ERROR "\tTraverse parameter NIHU_FMM_TRAVERSE must either be \"BFS\" or \"DFS\"")
endif()
list(APPEND NIHU_FMM_COMPILE_DEFINITIONS "-DNIHU_FMM_TRAVERSE_${NIHU_FMM_TRAVERSE}")

message(STATUS "\tFMM module configured successfully")
