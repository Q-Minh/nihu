CXX=g++
LFLAGS=-lgtest -lgtest_main -lpthread
CFLAGS=-std=c++11 -I../ -I/usr/local/include/eigen3

test_space.exe : test_space.o
	$(CXX) test_space.o $(LFLAGS) -o test_space.exe

test_space.o : test_space.cpp
	$(CXX) -c $(CFLAGS) test_space.cpp

