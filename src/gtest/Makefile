CXX=g++
TOOLBOXDIR=../../..
GTESTDIR=$(TOOLBOXDIR)/gtest
EIGENDIR=$(TOOLBOXDIR)/eigen
NIHUDIR=$(TOOLBOXDIR)/nihu
NIHULIBDIR=$(NIHUDIR)/src/library
LFLAGS=-lpthread
CFLAGS=-Wall -Wextra -std=c++11 -c -I../ -I../../../eigen -I$(GTESTDIR)/include -I$(GTESTDIR)

test_core.exe : test_space.o test_domain.o test_shape.o lib_domain.o lib_shape.o gtest_main.o gtest-all.o
	$(CXX) test_space.o test_domain.o test_shape.o lib_domain.o lib_shape.o gtest_main.o gtest-all.o $(LFLAGS) -o test_space.exe

gtest_main.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest_main.cc -o gtest_main.o

gtest-all.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest-all.cc -o gtest-all.o

test_space.o : test_space.cpp
	$(CXX) $(CFLAGS) test_space.cpp -o test_space.o

test_domain.o : test_domain.cpp
	$(CXX) $(CFLAGS) test_domain.cpp -o test_domain.o

test_shape.o : test_shape.cpp
	$(CXX) $(CFLAGS) test_shape.cpp -o test_shape.o

lib_domain.o :  $(NIHULIBDIR)/lib_domain.cpp
	$(CXX) $(CFLAGS) $(NIHULIBDIR)/lib_domain.cpp -o lib_domain.o

lib_shape.o :  $(NIHULIBDIR)/lib_shape.cpp
	$(CXX) $(CFLAGS) $(NIHULIBDIR)/lib_shape.cpp -o lib_shape.o

clean:
	rm *.o *~ *.exe

