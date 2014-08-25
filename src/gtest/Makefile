CXX=g++
TOOLBOXDIR=../../..
GTESTDIR=$(TOOLBOXDIR)/gtest
EIGENDIR=$(TOOLBOXDIR)/eigen
NIHUDIR=$(TOOLBOXDIR)/nihu
NIHULIBDIR=$(NIHUDIR)/src/library
LFLAGS=-lpthread
CFLAGS=-Wall -Wextra -std=c++11 -c -I../ -I../../../eigen -I$(GTESTDIR)/include -I$(GTESTDIR)

test_space.exe : test_space.o lib_domain.o gtest_main.o gtest-all.o
	$(CXX) test_space.o lib_domain.o gtest_main.o gtest-all.o $(LFLAGS) -o test_space.exe

gtest_main.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest_main.cc -o gtest_main.o

gtest-all.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest-all.cc -o gtest-all.o

#gtest_main.a: gtest_main.o gtest-all.o
#	ar rv gtest_main.a gtest-all.o gtest_main.o

test_space.o : test_space.cpp
	$(CXX) -c $(CFLAGS) test_space.cpp -o test_space.o

lib_domain.o :  $(NIHULIBDIR)/lib_domain.cpp
	$(CXX) -c $(CFLAGS) $(NIHULIBDIR)/lib_domain.cpp -o lib_domain.o

clean:
	rm *.o *~ *.exe

