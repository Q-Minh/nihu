SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -std=c++11" )

enable_testing()

add_subdirectory(tmp_unit)
add_subdirectory(core_unit)
