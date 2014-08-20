TOOLBOXDIR=../../..
GTESTDIR=$(TOOLBOXDIR)/gtest-1.7.0
EIGENDIR=$(TOOLBOXDIR)/eigen
NIHUDIR=$(TOOLBOXDIR)/nihu/src/nihu
NIHULIBDIR=$(TOOLBOXDIR)/nihu/src/library
CXX=g++
LFLAGS=-L$(GTESTDIR) -lgtest -lgtest_main -lpthread
CFLAGS=-std=c++11 -I../ -I../../../eigen -I$(GTESTDIR)/include

test_space.exe : test_space.o lib_domain.o
	$(CXX) test_space.o lib_domain.o $(LFLAGS) -o test_space.exe

test_space.o : test_space.cpp
	$(CXX) -c $(CFLAGS) test_space.cpp -o test_space.o

lib_domain.o :  $(NIHULIBDIR)/lib_domain.cpp
	$(CXX) -c $(CFLAGS) $(NIHULIBDIR)/lib_domain.cpp -o lib_domain.o
