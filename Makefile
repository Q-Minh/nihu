CC=g++
CFLAGS=-Wall -pedantic -O3 -I /usr/local/include/eigen3 -std=c++0x

element_test: src/bem/element_test.cpp src/bem/element.hpp src/bem/shapeset.hpp src/bem/element.hpp
	$(CC) $(CFLAGS) src/bem/element_test.cpp -o element_test

quad_test: src/bem/quadrature.cpp
	$(CC) $(CFLAGS) src/bem/quadrature.cpp -o quadrature_test

clean:
	rm element_test

