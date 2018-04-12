CXX:=g++
TOOLBOX_DIR:=../../../../..
NIHU_DIR:=$(TOOLBOX_DIR)/nihu/src
EIGEN_DIR:=$(TOOLBOX_DIR)/eigen
CFLAGS := -fopenmp -std=c++11 -Wall -c -I$(NIHU_DIR) -I$(EIGEN_DIR) -O4 -DITERATIVE
LFLAGS := -fopenmp

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
	
bem_const.exe: bem_const.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const.obj nihu.lib -o bem_const.exe

bem_gauss.exe: bem_gauss.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss.obj nihu.lib -o bem_gauss.exe

bem_const_bm.exe: bem_const_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_const_bm.obj nihu.lib -o bem_const_bm.exe

bem_gauss_bm.exe: bem_gauss_bm.obj nihu.lib
	$(CXX) $(LFLAGS) bem_gauss_bm.obj nihu.lib -o bem_gauss_bm.exe

all: bem_const.exe bem_gauss.exe bem_const_bm.exe bem_gauss_bm.exe

clean:
	rm main.exe *.obj *.lib
