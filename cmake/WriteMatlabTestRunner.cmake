# script for writing simple shell script that executes a matlab test

FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${run_file}" "${mat_root}/bin/matlab -nodisplay < ${m_test_file}")