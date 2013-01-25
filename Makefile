CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++0x

tests: quad_test element_test mesh_test field_test integral_test weighted_test rayleigh_test

quad_test: src/bem/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/bem/quad_test.cpp -o quad_test

element_test: src/bem/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/bem/element_test.cpp -o element_test

mesh_test: src/bem/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/bem/mesh_test.cpp -o mesh_test

field_test: src/bem/field_test.cpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/bem/field_test.cpp -o field_test

integral_test: src/bem/integral_test.cpp src/bem/kernel.hpp src/bem/integral.hpp
	$(CC) $(CFLAGS) src/bem/integral_test.cpp -o integral_test

weighted_test: src/bem/weighted_integral_test.cpp src/bem/weighted_integral.hpp
	$(CC) $(CFLAGS) src/bem/weighted_integral_test.cpp -o weighted_test

rayleigh_test: src/bem/rayleigh_test.cpp src/bem/rayleigh.hpp src/bem/weighted_integral.hpp  src/bem/function_space.hpp
	$(CC) $(CFLAGS) src/bem/rayleigh_test.cpp -o rayleigh_test

bem_test: src/bem/bem_test.cpp src/bem/bem.hpp src/bem/weighted_integral.hpp  src/bem/function_space.hpp
	$(CC) $(CFLAGS) src/bem/bem_test.cpp -o bem_test

clean:
	rm *_test

