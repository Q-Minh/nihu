include make.inc

CXX:=g++
CFLAGS := -fopenmp -std=c++11 -Wall -c -I$(NIHU_DIR) -I$(EIGEN_DIR) -O4 -DITERATIVE -DUSE_GMRES -DNUM_PROCESSORS=1
LFLAGS := -fopenmp

%.obj: $(NIHU_DIR)/library/%.cpp
	$(CXX) $(CFLAGS) $< -o $@

nihu.lib: lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj
	ar rcs nihu.lib lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj

bem_const.obj: bem.cpp
	$(CXX) $(CFLAGS) bem.cpp -o $@
	
bem_const_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) bem.cpp -DBURTON_MILLER -o $@
	
bem_gauss.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS bem.cpp -o $@
	
bem_gauss_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS -DBURTON_MILLER bem.cpp -o $@
	
bem_const.exe: bem_const.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const.obj nihu.lib -o $@

bem_gauss.exe: bem_gauss.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss.obj nihu.lib -o $@

bem_const_bm.exe: bem_const_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const_bm.obj nihu.lib -o $@

bem_gauss_bm.exe: bem_gauss_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss_bm.obj nihu.lib -o $@

all: bem_const.exe bem_gauss.exe bem_const_bm.exe bem_gauss_bm.exe

clean:
	rm main.exe *.obj *.lib
