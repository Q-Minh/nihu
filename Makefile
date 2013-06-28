CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++11

space_test: src/test/space_test.cpp src/bem/space.hpp
	$(CC) $(CFLAGS) src/test/space_test.cpp -o space_test

domain_test: src/test/domain_test.cpp src/bem/domain.hpp src/bem/space.hpp 
	$(CC) $(CFLAGS) src/test/domain_test.cpp -o domain_test

shape_test: src/test/shape_test.cpp src/bem/shapeset.hpp src/bem/domain.hpp
	$(CC) $(CFLAGS) src/test/shape_test.cpp -o shape_test

element_test: src/test/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp
	$(CC) $(CFLAGS) src/test/element_test.cpp -o element_test

mesh_test: src/test/mesh_test.cpp src/bem/mesh.hpp
	$(CC) $(CFLAGS) src/test/mesh_test.cpp -o mesh_test

field_test: src/test/field_test.cpp src/bem/field.hpp
	$(CC) $(CFLAGS) src/test/field_test.cpp -o field_test

quad_test: src/test/quad_test.cpp src/bem/quadrature.hpp
	$(CC) $(CFLAGS) src/test/quad_test.cpp -o quad_test

weighted_residual_test: src/test/weighted_residual_test.cpp
	$(CC) $(CFLAGS) src/test/weighted_residual_test.cpp  -o weighted_residual_test

weighted_residual_test_q: src/test/weighted_residual_test_q.cpp
	$(CC) $(CFLAGS) src/test/weighted_residual_test_q.cpp  -o weighted_residual_test_q

field_accelerator_test: src/test/field_accelerator_test.cpp src/bem/field_type_accelerator.hpp
	$(CC) $(CFLAGS) src/test/field_accelerator_test.cpp -o field_accelerator_test

couple_test: src/test/couple_test.cpp src/bem/couple.hpp
	$(CC) $(CFLAGS) src/test/couple_test.cpp -o couple_test

singulartity_check_test: src/test/singularity_check_test.cpp src/bem/element.hpp src/bem/singularity_check.hpp
	$(CC) $(CFLAGS) src/test/singularity_check_test.cpp -o singularity_check_test

gaussian_test: examples/gaussian_test.cpp
	$(CC) $(CFLAGS) examples/gaussian_test.cpp -o gaussian_test

clean:
	rm *_test

