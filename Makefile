CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++0x

tests: quad_test element_test mesh_test field_test

quad_test: src/bem/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/bem/quad_test.cpp -o quad_test

element_test: src/bem/element_test.cpp src/bem/elem_descriptor.hpp src/bem/element.hpp src/bem/shapeset.hpp src/bem/element.hpp
	$(CC) $(CFLAGS) src/bem/element_test.cpp -o element_test

mesh_test: src/bem/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/bem/mesh_test.cpp -o mesh_test

field_test: src/bem/field_test.cpp src/bem/field.hpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/bem/field_test.cpp -o field_test

clean:
	rm *test

