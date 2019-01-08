include make.inc

CXX:=g++
CFLAGS := -fopenmp -std=c++11 -Wall -c -I$(NIHU_DIR) -I$(EIGEN_DIR) -O4 -DITERATIVE
LFLAGS := -fopenmp

TARGETS:= bem_const bem_gauss bem_const_bm bem_gauss_bm

all: $(TARGETS)

%.obj: $(NIHU_DIR)/library/%.cpp
	$(CXX) $(CFLAGS) $< -o $@

nihu.lib: lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj
	ar rcs nihu.lib lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj

bem_const.obj: bem.cpp
	$(CXX) $(CFLAGS) bem.cpp -o bem_const.obj
	
bem_const_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) bem.cpp -DBM -o bem_const_bm.obj
	
bem_gauss.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS bem.cpp -o bem_gauss.obj
	
bem_gauss_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS -DBM bem.cpp -o bem_gauss_bm.obj
	
bem_const: bem_const.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const.obj nihu.lib -o $@

bem_gauss: bem_gauss.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss.obj nihu.lib -o $@

bem_const_bm: bem_const_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const_bm.obj nihu.lib -o $@

bem_gauss_bm: bem_gauss_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss_bm.obj nihu.lib -o $@

clean:
	rm -f $(TARGETS) *.obj *.lib
