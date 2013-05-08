CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++11

domain.o:  src/bem/domain.cpp 
	$(CC) -c $(CFLAGS) src/bem/domain.cpp -o domain.o

shapeset.o: src/bem/shapeset.cpp
	$(CC) -c $(CFLAGS) src/bem/shapeset.cpp -o shapeset.o


quad_test: src/test/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/test/quad_test.cpp -o quad_test

element_test: src/test/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/test/element_test.cpp -o element_test

mesh_test: src/test/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/test/mesh_test.cpp -o mesh_test

field_test: src/test/field_test.cpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/test/field_test.cpp -o field_test

weighted_residual_test: domain.o shapeset.o src/test/weighted_residual_test.cpp
	$(CC) -c $(CFLAGS) src/test/weighted_residual_test.cpp -o weighted_residual_test.o
	$(CC) weighted_residual_test.o domain.o shapeset.o -o weighted_residual_test

shape_test: src/test/shape_test.cpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/test/shape_test.cpp -o shape_test

field_accelerator_test: src/test/field_accelerator_test.cpp src/bem/field_type_accelerator.hpp
	$(CC) $(CFLAGS) src/test/field_accelerator_test.cpp -o field_accelerator_test

couple_test: src/test/couple_test.cpp src/bem/couple.hpp
	$(CC) $(CFLAGS) src/test/couple_test.cpp -o couple_test

singulartity_check_test: src/test/singularity_check_test.cpp src/bem/element.hpp src/bem/singularity_check.hpp
	$(CC) $(CFLAGS) src/test/singularity_check_test.cpp -o singularity_check_test

clean:
	rm *_test

