@ECHO OFF
SET CFLAGS=-O3 -Wall -std=c++11 -I../../../Eigen
@ECHO ON
g++ %CFLAGS% space_test.cpp -o space_test.exe
g++ %CFLAGS% domain_test.cpp -o domain_test.exe
g++ %CFLAGS% shape_test.cpp -o shape_test.exe
g++ %CFLAGS% element_test.cpp -o element_test.exe
