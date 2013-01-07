#include "mesh.hpp"

#include <algorithm>
#include <iostream>


int main(void)
{
	Mesh<3> mesh;
	
	mesh.add_elem(QuadElem());
	mesh.add_elem(TriaElem());
	
	std::for_each(
		mesh.begin<QuadElem>(),
		mesh.end<QuadElem>(),
		[] (QuadElem &e) { std::cout << e[0] << std::endl; }
	);

	std::for_each(
		mesh.begin<QuadElem>(),
		mesh.end<QuadElem>(),
		[] (QuadElem &e) { std::cout << e[0] << std::endl; }
	);

	return 0;
}

