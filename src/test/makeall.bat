@ECHO OFF
SET CFLAGS=-O3 -Wall -pedantic -std=c++11 -I../../../Eigen -I..
@ECHO ON
g++ %CFLAGS% core_unit\space_test.cpp -o space_test.exe
g++ %CFLAGS% core_unit\domain_test.cpp -o domain_test.exe
g++ %CFLAGS% core_unit\shape_test.cpp -o shape_test.exe
g++ %CFLAGS% core_unit\element_test.cpp -o element_test.exe
g++ %CFLAGS% core_unit\field_test.cpp -o field_test.exe
