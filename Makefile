CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++0x

tests: quad_test element_test mesh_test field_test weighted_residual_test shape_test field_accelerator_test couple_test

quad_test: src/test/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/test/quad_test.cpp -o quad_test

element_test: src/test/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/test/element_test.cpp -o element_test

mesh_test: src/test/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/test/mesh_test.cpp -o mesh_test

field_test: src/test/field_test.cpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/test/field_test.cpp -o field_test

weighted_residual_test: src/test/weighted_residual_test.cpp src
	$(CC) $(CFLAGS) src/test/weighted_residual_test.cpp -o weighted_residual_test

shape_test: src/test/shape_test.cpp src/bem/shapeset.hpp src
	$(CC) $(CFLAGS) src/test/shape_test.cpp -o shape_test

field_accelerator_test: src/test/field_accelerator_test.cpp src/bem/field_type_accelerator.hpp src
	$(CC) $(CFLAGS) src/test/field_accelerator_test.cpp -o field_accelerator_test

couple_test: src/test/couple_test.cpp src/bem/couple.hpp src
	$(CC) $(CFLAGS) src/test/couple_test.cpp -o couple_test

clean:
	rm *_test

