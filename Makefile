CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++0x

tests: quad_test element_test mesh_test field_test weighted_test rayleigh_test bem_test

quad_test: src/test/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/test/quad_test.cpp -o quad_test

element_test: src/test/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/test/element_test.cpp -o element_test

mesh_test: src/test/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/test/mesh_test.cpp -o mesh_test

field_test: src/test/field_test.cpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/test/field_test.cpp -o field_test

rayleigh_test: src/test/rayleigh_test.cpp src/bem/bem.hpp src/bem/weighted_integral.hpp  src/bem/function_space.hpp
	$(CC) $(CFLAGS) src/test/rayleigh_test.cpp -o rayleigh_test

bem_test: src/test/bem_test.cpp src/bem/bem.hpp src/bem/bem.hpp src/bem/weighted_integral.hpp  src/bem/function_space.hpp
	$(CC) $(CFLAGS) src/test/bem_test.cpp -o bem_test

accelerator_test: src/test/accelerator_test.cpp src
	$(CC) $(CFLAGS) src/test/accelerator_test.cpp -o accelerator_test

double_integral_test: src/test/surface_test.cpp src
	$(CC) $(CFLAGS) src/test/surface_test.cpp -o double_integral_test

weighted_residual_test: src/test/weighted_residual_test.cpp src
	$(CC) $(CFLAGS) src/test/weighted_residual_test.cpp -o weighted_residual_test

clean:
	rm *_test

