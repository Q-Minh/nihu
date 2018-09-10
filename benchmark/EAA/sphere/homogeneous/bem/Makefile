CXX:=g++
TOOLBOX_DIR:=d:/toolbox
NIHU_DIR:=$(TOOLBOX_DIR)/nihu/src
EIGEN_DIR:=$(TOOLBOX_DIR)/eigen
MATLAB_DIR:="c:/Program Files/MATLAB/R2016a"
CFLAGS := -fopenmp -std=c++11 -Wall -c -I$(NIHU_DIR) -I$(EIGEN_DIR) -I$(MATLAB_DIR)/extern/include -O4 -DITERATIVE
LFLAGS := -fopenmp

%.obj: $(NIHU_DIR)/library/%.cpp
	$(CXX) $(CFLAGS) $< -o $@

nihu.lib: lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj
	ar rcs nihu.lib lib_domain.obj lib_shape.obj lib_element.obj lib_tmp.obj lib_singular_galerkin_quadrature.obj

bem_const.obj: bem.cpp
	$(CXX) $(CFLAGS) -DKHIE bem.cpp -o bem_const.obj
	
bem_const_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) -DBM bem.cpp -o bem_const_bm.obj
	
bem_const_hs.obj: bem.cpp
	$(CXX) $(CFLAGS) -DHSIE bem.cpp -o bem_const_hs.obj
	
bem_gauss.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS -DKHIE bem.cpp -o bem_gauss.obj
	
bem_gauss_bm.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS -DBM bem.cpp -o bem_gauss_bm.obj
	
bem_gauss_hs.obj: bem.cpp
	$(CXX) $(CFLAGS) -DGAUSS -DHSIE bem.cpp -o bem_gauss_hs.obj
	
	
%.mex.obj: %.mex.cpp
	$(CXX) $(CFLAGS) $< -o $@
	

%.exe: %.obj nihu.lib
	$(CXX) $(LFLAGS) $^ -o $@
	
	
%.mexw64: %.mex.obj nihu.lib
	$(CXX) -shared -m64 -L$(MATLAB_DIR)/bin/win64 -L$(MATLAB_DIR)/extern/lib/win64/microsoft -L$(FFTW_DIR) -L. $< nihu.lib -lmex -lmx -lmat -o $@
	
	

all: bem_const.exe bem_gauss.exe bem_const_bm.exe bem_gauss_bm.exe bem_const_hs.exe bem_gauss_hs.exe

mex: bem_matrices.mexw64


clean:
	del *.exe *.obj *.lib
