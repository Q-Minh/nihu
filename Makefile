CC=g++
CFLAGS=-Wall -pedantic -O3 -I/usr/local/include/eigen3 -std=c++0x

main: src/main.cpp
	$(CC) $(CFLAGS) src/main.cpp -o main
main_seq: src/main_seq.cpp
	$(CC) $(CFLAGS) src/main_seq.cpp -o main_seq

mesh: src/bem/mesh.cpp src/bem/mesh.hpp src/bem/node.hpp src/bem/element.hpp
	$(CC) $(CFLAGS) src/bem/mesh.cpp -o mesh

clean:
	rm main_seq main mesh
