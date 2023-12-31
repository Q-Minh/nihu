# script for writing simple shell script that executes a matlab test

if(WIN32)
	set(matlab_exec "matlab.exe")
	set(run_file "${run_file}.bat")
else(WIN32)
	set(matlab_exec "matlab")
endif(WIN32)

SET(write_str "${mat_root}/bin/${matlab_exec} -nodisplay < ${m_test_file}")

if(WIN32)
	STRING(REGEX REPLACE "/" "\\\\" write_str ${write_str})
endif(WIN32)

FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${run_file}" "${write_str}")