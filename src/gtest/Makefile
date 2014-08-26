TOOLBOXDIR := ../../..
GTESTDIR := $(TOOLBOXDIR)/gtest
EIGENDIR := $(TOOLBOXDIR)/eigen
NIHUSRCDIR := $(TOOLBOXDIR)/nihu/src
NIHULIBDIR := $(NIHUSRCDIR)/library
TEST_CPP_FILES := $(wildcard *.cpp)
TEST_OBJ_FILES := $(TEST_CPP_FILES:.cpp=.o)
LIB_CPP_FILES := $(wildcard $(NIHULIBDIR)/*.cpp)
LIB_OBJ_FILES := $(LIB_CPP_FILES:.cpp=.o)
GTEST_CPP_FILES := $(GTESTDIR)/src/gtest_main.cc $(GTESTDIR)/src/gtest-all.cc
GTEST_OBJ_FILES := $(GTEST_CPP_FILES:.cc=.o)
LFLAGS := -lpthread
CFLAGS := -Wall -Wextra -std=c++11 -c -isystem $(GTESTDIR)/include -I$(NIHUSRCDIR) -I$(EIGENDIR) -I$(GTESTDIR)

test_main.exe : $(TEST_OBJ_FILES) $(LIB_OBJ_FILES) $(GTEST_OBJ_FILES)
	$(CXX) -o $@ $^ $(LFLAGS)

%.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ $<

%.o : %.cc
	$(CXX) $(CFLAGS) -o $@ $<

clean:
	rm $(TEST_OBJ_FILES) $(GTEST_OBJ_FILES) $(LIB_OBJ_FILES)

