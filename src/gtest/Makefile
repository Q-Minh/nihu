CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
TOOLBOXDIR := ../../..
GTESTDIR := $(TOOLBOXDIR)/gtest
EIGENDIR := $(TOOLBOXDIR)/eigen
NIHUDIR := $(TOOLBOXDIR)/nihu
NIHULIBDIR := $(NIHUDIR)/src/library
LFLAGS := -lpthread
CFLAGS := -Wall -Wextra -std=c++11 -c -I../ -I../../../eigen -I$(GTESTDIR)/include -I$(GTESTDIR)

test_main.exe : $(OBJ_FILES) lib_domain.o lib_shape.o gtest_main.o gtest-all.o
	$(CXX) -o $@ $^ $(LFLAGS)

gtest_main.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest_main.cc -o gtest_main.o

gtest-all.o:
	$(CXX) $(CFLAGS) -isystem $(GTESTDIR)/include -I$(GTESTDIR) $(GTESTDIR)/src/gtest-all.cc -o gtest-all.o

%.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ $<

lib_domain.o :  $(NIHULIBDIR)/lib_domain.cpp
	$(CXX) $(CFLAGS) $(NIHULIBDIR)/lib_domain.cpp -o lib_domain.o

lib_shape.o :  $(NIHULIBDIR)/lib_shape.cpp
	$(CXX) $(CFLAGS) $(NIHULIBDIR)/lib_shape.cpp -o lib_shape.o

clean:
	rm *.o *~ *.exe

